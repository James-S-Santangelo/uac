# Predicting the strength of urban-rural clines in a 
# Mendelian polymorphism along a latitudinal gradient 
#
# Authors: James S. Santangelo, Ken A. Thompson, Beata Cohan
# Jibran Syed, Rob W. Ness, Marc T. J. Johnson
#
#
# Generates tables and figures for manuscript

###########################
#### TABLES: MAIN TEXT ####
###########################

## TABLE 1 

# Linear cline models only. Part of table 1
linearClines <- read_csv("analysis/clineModelOutput.csv") %>% 
  select(City, betaLinOnly, betaLin_AcHWE, betaLin_LiHWE)

# Extract number of populations and plants per city. Merge with linear clines model output above
table1 <- datPlants %>%
  group_by(City) %>%
  summarise(numPops = n_distinct(Population),
            numPlants = n()) %>%
  merge(., linearClines, by = "City")

# Write table 1 to disk
write_csv(table1, "analysis/tables/main-text/Table1_cityClineSummary.csv")

###############################
#### TABLES: SUPPLEMENTARY ####
###############################

## TABLE S1

weather_data <- read_csv("data-clean/DailyNormals_AllCities_Filtered.csv")

# Summarize number of observations for each city
weather_summ <- weather_data %>% 
  group_by(City, STATION_NAME) %>%
  summarise(min_year = min(Year),
            max_year = max(Year),
            count = n()) %>%
  mutate(STATION_NAME = str_to_title(STATION_NAME))

# Write summarized data to disk
write_csv(weather_summ, "analysis/tables/supplemental/TableS1_DailyNormals_Summary.csv")

## TABLE S2

# Write pairwise correlation matrix to disk
write_csv(corr_mat, "analysis/tables/supplemental/TableS2_weatherCorrMat.csv")

## TABLE S3

# Get the best fit cline model for each city and locus
# Possible values for the response variables
possibleResponses <- c("freqHCN", "AcHWE", "LiHWE")

tbl_out <- c() # vector to hold resulting dataframes

# Loop over response variables
for(res in possibleResponses){

  tbl_out[[res]] <- datPops %>% 
    split(.$City) %>%  # Separate list for each city
    map_dfr(., ~tidy(getBestFitClineModelOrder(res, .)), .id = "City") %>% 
    mutate(response = res) %>% 
    spread(key = response, value = x)
  
}

# Left join all dataframes with model orders
allModelOutputs <- reduce(tbl_out, left_join, by = "City")

# Write model order table to disk
write_csv(allModelOutputs, "analysis/tables/supplemental/TableS3_clineOrderData.csv")

## TABLE S4

# Write HCN dredge models to disk
write_csv(HCN_dredge_models, path = "analysis/tables/supplemental/TableS4_HCN_dredge_output.csv", 
          col_names = TRUE)

## TABLE S5

# Write full model coefficients from dredge output to disk
tableS5 <- as.data.frame(rbind(
  c("Full"),
  summary(HCN_modAvg)$coefmat.full)) %>%
  rownames_to_column()

write_csv(tableS5, "analysis/tables/supplemental/TableS5_freqHCN_FullModelAvg.csv")

############################
#### FIGURES: MAIN TEXT ####
############################

#Theme used to plot figures throughout script.
ng1 = theme(
  aspect.ratio = 0.7,
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(color = "black", size = 1),
  axis.line.y = element_line(color = "black", size = 1),
  axis.ticks = element_line(color = "black"),
  axis.text = element_text(color = "black", size = 15),
  axis.title = element_text(color = "black", size = 1),
  axis.title.y = element_text(vjust = 2, face = "bold", size = 18),
  axis.title.x = element_text(vjust = 0.1, face = "bold", size = 18),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  legend.position = "none",
  legend.direction = "vertical",
  legend.text = element_text(size = 13),
  legend.key = element_rect(fill = "white"),
  legend.title = element_text(size = 15, face = "bold"),
  legend.key.size = unit(1.0, "cm")
)

## FIGURE 1

# Figure 1 was generated in QGIS (v. 3.2.3) and is not included here

## FIGURE 2 

# Plot of HCN frequency against distance for each city
# Solid line if significant. Thick black line is regression across all cities.
HCN_by_city <- datPops %>%
  group_by(City) %>%
  do(mod = lm(freqHCN ~ std_distance, data = .)) %>% # Linear model for each city
  tidy(., mod) %>% # Clean up coefficients
  mutate(intercept = estimate[term == "(Intercept)"]) %>%
  filter(term == "std_distance") %>% # Only need beta coefficient
  mutate(predicted = intercept + (estimate * 1.0)) %>% # Calculate predicted value along regression line for each city
  select(City, p.value, predicted) %>%
  merge(., datPops, by = "City", all.y = TRUE) %>% # Merge predicted values for each city to pop-mean dataframe
  mutate(significant = ifelse(p.value < 0.05, "Yes", "No")) %>%
  ggplot(., aes(x = std_distance, y = freqHCN)) +  
  # Regression line for each city
  geom_line(stat = "smooth", method="lm", 
            aes(linetype = significant,
                group = City), 
            alpha = 0.7, size = 1) +
  # Overall regression line
  geom_line(stat = "smooth", method="lm", colour = "black", size = 2.5) +
  # scale_colour_manual(values = colors) +
  xlab("Standardized distance") + ylab("Frequency of HCN") + 
  scale_linetype_manual(values=c("dashed", "solid")) +
  scale_y_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
  scale_x_continuous(breaks = seq(from = 0, to = 1.1, by = 0.25)) +
  coord_cartesian(ylim = c(0.1, 1.025), xlim = c(0, 1), clip = 'off') +
  geom_text(aes(label = City, x = 1.005, y = predicted), hjust = 0) + 
  ng1 +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.key.height = unit(0.5, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 2)))
# geom_dl(aes(label = City), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.8))
HCN_by_city

# Save file to disk
ggsave(filename = "analysis/figures/main-text/Figure2_HCN-by-distance.pdf", 
       plot = HCN_by_city, device = 'pdf', units = 'in',
       width = 12, height = 8, dpi = 600)



## FIGURE 3 ##

# Add abbreviations
citySummaryData <- citySummaryData %>%
  mutate(abbr = case_when(City == "Jacksonville" ~ "Jax",
                          City == "Tampa" ~ "Tpa",
                          City == "Atlanta" ~ "Atl",
                          City == "Norfolk" ~ "Nor",
                          City == "Charlotte" ~ "Clt",
                          City == "Toronto" ~ "Tor",
                          City == "Montreal" ~ "Mtl",
                          City == "Detroit" ~ "Det",
                          City == "Washington D.C." ~ "DC",
                          City == "Cleveland" ~ "Clv",
                          City == "NewYork" ~ "NY",
                          City == "Pittsburgh" ~ "Pgh",
                          City == "Boston" ~ "Bos",
                          City == "Baltimore" ~ "Blt",
                          City == "Cincinnati" ~ "Cin",
                          TRUE ~ "Phl"))


# Figure 3a. HCN against # days < 0 with no snow
HCN_by_DaysNeg <- citySummaryData %>%
  ggplot(., aes(x = daysNegNoSnow, y = freqHCN)) +
  # geom_point(size = 2.5) +
  geom_smooth(method = "lm", size = 1.5, colour = "black", 
              se = FALSE) +
  scale_x_continuous(breaks = seq(from = 0, to = 35, by = 5)) +
  xlab("# days < 0°C with no snow") + ylab("Mean HCN frequency") +
  geom_text(aes(label = abbr), vjust = 0, hjust = 0) + 
  ng1
HCN_by_DaysNeg

ggsave(filename = "analysis/figures/main-text/Figure3a_HCN-by-NumDaysNegNoSnow.pdf", 
       plot = HCN_by_DaysNeg, device = 'pdf', units = 'in',
       width = 5, height = 5, dpi = 600)

# Figure 3b. HCN against PC1
HCN_by_PC1 <- citySummaryData %>%
  ggplot(., aes(x = PC1_HCN, y = freqHCN)) +
  # geom_point(size = 2.5) +
  geom_smooth(method = "lm", size = 1.5, colour = "black", 
              se = FALSE) +
  # scale_x_continuous(breaks = seq(from = 0, to = 35, by = 5)) +
  xlab("PC1 (90.2%)") + ylab("Mean HCN frequency") +
  geom_text(aes(label = abbr), vjust = 0, hjust = 0) + 
  ng1
HCN_by_PC1

ggsave(filename = "analysis/figures/main-text/Figure3b_HCN-by-PC1.pdf", 
       plot = HCN_by_PC1, device = 'pdf', units = 'in',
       width = 5, height = 5, dpi = 600)


## FIGURE 4 

citySummaryDataForAnalysis <- citySummaryDataForAnalysis %>% 
  left_join(., citySummaryData %>% select(City, abbr), by = "City")

# Figure 4. Slope against PC1
Slope_by_PC1Lin <- citySummaryDataForAnalysis %>%
  ggplot(., aes(x = PC1_SlopeLin, y = betaLinOnly)) +
  # geom_point(size = 2.5) +
  geom_smooth(method = "lm", size = 1.5, colour = "black", 
              se = FALSE) +
  # scale_x_continuous(breaks = seq(from = 0, to = 35, by = 5)) +
  xlab("PC1 (92.8%)") + ylab("Slope of HCN cline") +
  geom_text(aes(label = abbr), vjust = 0, hjust = 0) + 
  ng1
Slope_by_PC1Lin

ggsave(filename = "analysis/figures/main-text/Figure4_Slope-by-PC1.pdf", 
       plot = Slope_by_PC1Lin, device = 'pdf', units = 'in',
       width = 5, height = 5, dpi = 600)

################################
#### FIGURES: SUPPLEMENTARY ####
################################

## FIGURE S1

# Figure S1 is generated with scripts/powerAnalysis.R

## FIGURE S2

#Plot of HCN frequencies with latitude
plotHCN_by_lat <- ggplot(citySummaryData, aes(x = Latitude, y = freqHCN)) +
  geom_point(colour = "black", size = 3.5) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", size = 2) + 
  ylab("Frequency of HCN") + xlab("Latitude") +
  ng1
plotHCN_by_lat

ggsave(filename = "analysis/figures/supplemental/FigureS2_HCN_by_Lat.pdf", 
       plot = plotHCN_by_lat, device = "pdf", 
       width = 5, height = 5, dpi = 300)

## FIGURES S5 - S20 (INDIVIDUAL CLINE BIPLOTS)

# Create list with city dataframes as elements
city_df_list <- datPops %>% split(.$City)

# Create biplot for all cities with HCN as response
outpath <- "analysis/figures/individualCline_biplots/HCN/"
purrr::walk(city_df_list, clineBiplot, 
            response_var = "freqHCN", 
            outpath = outpath,
            model_order_df = clineModelOrder)

# Create biplot for all cities with Ac as response
outpath <- "analysis/figures/individualCline_biplots/Ac/"
purrr::walk(city_df_list, clineBiplot, 
            response_var = "AcHWE", 
            outpath = outpath,
            model_order_df = clineModelOrder)

# Create biplot for all cities with Li as response
outpath <- "analysis/figures/individualCline_biplots/Li/"
purrr::walk(city_df_list, clineBiplot, 
            response_var = "LiHWE", 
            outpath = outpath,
            model_order_df = clineModelOrder)

## FIGURE SXX ##

# HCN by city logistic
HCN_by_cityLog <- datPlants %>%
  group_by(City) %>%
  do(mod = glm(HCN_Result ~ std_distance, family = "binomial", data = .)) %>%
  tidy(., mod) %>%
  mutate(intercept = estimate[term == "(Intercept)"]) %>%
  filter(term == "std_distance") %>%
  mutate(predicted = intercept + (estimate * 1.0),
         odds = exp(predicted),
         prob = odds / (1 + odds)) %>%
  select(City, p.value, predicted, odds, prob) %>%
  merge(., datPlants, by = "City", all.y = TRUE) %>%
  mutate(significant = ifelse(p.value < 0.05, "Yes", "No")) %>%
  ggplot(., aes(x = std_distance, y = HCN_Result)) +  
  geom_line(stat = "smooth", method="glm", method.args = list(family = "binomial"),
            aes(linetype = significant,
                group = City), 
            alpha = 0.7, size = 1) +
  geom_line(stat = "smooth", method="glm", method.args = list(family = "binomial"),
            colour = "black", size = 2.5) +
  # scale_colour_manual(values = colors) +
  xlab("Standardized distance") + ylab("Probability of being HCN+") + 
  scale_linetype_manual(values=c("dashed", "solid")) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25)) +
  scale_x_continuous(breaks = seq(from = 0, to = 1.1, by = 0.25)) +
  coord_cartesian(ylim = c(0, 1.025), xlim = c(0, 1), clip = 'off') +
  geom_text(aes(label = City, x = 1.005, y = prob), hjust = 0) + 
  ng1 +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.key.height = unit(0.5, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 2)))
# geom_dl(aes(label = City), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.8))
HCN_by_cityLog

## FIGURE SXX ##

# 4 supplementary multi-paneled figures with sigmoidal clines
# for each city

# Figure SX
ATL_LogReg <- plotLogReg(datPlants, "Atlanta", tag = "(a)")
BTL_LogReg <- plotLogReg(datPlants, "Baltimore", tag = "(b)")
BOS_LogReg <- plotLogReg(datPlants, "Boston", tag = "(c)")
CLT_LogReg <- plotLogReg(datPlants, "Charlotte", tag = "(d)")

logRegs1 <- ATL_LogReg + BTL_LogReg + BOS_LogReg + CLT_LogReg + plot_layout(ncol = 2)

ggsave(filename = "analysis/figures/supplemental/figureSX_logRegs_ATL-CLT.pdf", 
       plot = logRegs1, device = "pdf", 
       width = 8, height = 8, units = "in", dpi = 300) 

# Figure SX
CIN_LogReg <- plotLogReg(datPlants, "Cincinnati", tag = "(a)")
CLE_LogReg <- plotLogReg(datPlants, "Cleveland", tag = "(b)")
DET_LogReg <- plotLogReg(datPlants, "Detroit", tag = "(c)")
JAX_LogReg <- plotLogReg(datPlants, "Jacksonville", tag = "(d)")

logRegs2 <- CIN_LogReg + CLE_LogReg + DET_LogReg + JAX_LogReg + plot_layout(ncol = 2)

ggsave(filename = "analysis/figures/supplemental/figureSX_logRegs_CIN-JAX.pdf", 
       plot = logRegs2, device = "pdf", 
       width = 8, height = 8, units = "in", dpi = 300) 

# Figure SX
MTL_LogReg <- plotLogReg(datPlants, "Montreal", tag = "(a)")
NY_LogReg <- plotLogReg(datPlants, "NewYork", tag = "(b)")
NOR_LogReg <- plotLogReg(datPlants, "Norfolk", tag = "(c)")
PHL_LogReg <- plotLogReg(datPlants, "Philadelphia", tag = "(d)")

logRegs3 <- MTL_LogReg + NY_LogReg + NOR_LogReg + PHL_LogReg + plot_layout(ncol = 2)

ggsave(filename = "analysis/figures/supplemental/figureSX_logRegs_MTL-PHL.pdf", 
       plot = logRegs3, device = "pdf", 
       width = 8, height = 8, units = "in", dpi = 300) 

# Figure SX
PIT_LogReg <- plotLogReg(datPlants, "Pittsburgh", tag = "(a)")
# TMP_LogReg <- plotLogReg(datPlants, "Tampa", tag = "(n)")
TOR_LogReg <- plotLogReg(datPlants, "Toronto", tag = "(b)")
WDC_LogReg <- plotLogReg(datPlants, "Washington D.C.", tag = "(c)")

logRegs4 <- PIT_LogReg + TOR_LogReg + WDC_LogReg + plot_layout(ncol = 2)

ggsave(filename = "analysis/figures/supplemental/figureSX_logRegs_PIT-WDC.pdf", 
       plot = logRegs4, device = "pdf", 
       width = 8, height = 8, units = "in", dpi = 300) 

## FIGURE SXX ##

plotBosSlopePlants <- SlopeBosPlants %>% 
  select_if(is.numeric) %>% 
  purrr::map_dfr(~data.frame(prob = calcProb(., obs_BosSlope, nreps))) %>% 
  mutate(num_plants = 20:10) %>% 
  ggplot(., aes(x = num_plants, y = prob)) +
  geom_point(size = 2, colour = "black") +
  coord_cartesian(ylim = c(0.495, 0.505)) +
  scale_x_reverse(breaks = seq(from = 10, to = 20, by = 2)) +
  # scale_y_continuous(breaks = seq(from = 0.9, to = 0.98, by = 0.02)) +
  geom_vline(xintercept = 15, linetype = "dashed") +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  xlab("Number of plants") + ylab("P(slope \u2265 observed Boston slope)") +
  ng1
plotBosSlopePlants

# Save plot to disk
ggsave("analysis/figures/supplemental/samplingPlants_slopeBoston.pdf",
       plot = plotBosSlopePlants, dpi = 300, width = 6, height = 6, units = "in",
       device = cairo_pdf)