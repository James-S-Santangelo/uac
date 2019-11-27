# Predicting the strength of urban-rural clines in a 
# Mendelian polymorphism along a latitudinal gradient 
#
# Authors: James S. Santangelo, Ken A. Thompson, Beata Cohan
# Jibran Syed, Rob W. Ness, Marc T. J. Johnson
#
#
# Script to carry out all analyses reported in manuscript


### SETUP ###

# Change default contrasts to enable type III SS
options(contrasts = c("contr.sum", "contr.poly"))
# options(contrasts = c("contr.treatment", "contr.poly")) # Default

# Clear environement, if necessary
rm(list=ls())
.rs.restartR()

# Load required packages
library(tidyverse)
library(broom)
library(Hmisc)
library(factoextra)
library(vegan)
library(MuMIn)
library(car)
library(patchwork)
source('scripts/functions.R')

ATL_LogReg <- plotLogReg(allPlants, "Atlanta", tag = "(a)")
BTL_LogReg <- plotLogReg(allPlants, "Baltimore", tag = "(b)")
BOS_LogReg <- plotLogReg(allPlants, "Boston", tag = "(c)")
CLT_LogReg <- plotLogReg(allPlants, "Charlotte", tag = "(d)")


CIN_LogReg <- plotLogReg(allPlants, "Cincinnati", tag = "(a)")
CLE_LogReg <- plotLogReg(allPlants, "Cleveland", tag = "(b)")
DET_LogReg <- plotLogReg(allPlants, "Detroit", tag = "(c)")
JAX_LogReg <- plotLogReg(allPlants, "Jacksonville", tag = "(d)")


MTL_LogReg <- plotLogReg(allPlants, "Montreal", tag = "(a)")
NY_LogReg <- plotLogReg(allPlants, "NewYork", tag = "(b)")
NOR_LogReg <- plotLogReg(allPlants, "Norfolk", tag = "(c)")
PHL_LogReg <- plotLogReg(allPlants, "Philadelphia", tag = "(d)")


PIT_LogReg <- plotLogReg(allPlants, "Pittsburgh", tag = "(a)")
# TMP_LogReg <- plotLogReg(allPlants, "Tampa", tag = "(n)")
TOR_LogReg <- plotLogReg(allPlants, "Toronto", tag = "(b)")
WDC_LogReg <- plotLogReg(allPlants, "Washington D.C.", tag = "(c)")

logRegs1 <- ATL_LogReg + BTL_LogReg + BOS_LogReg + CLT_LogReg + plot_layout(ncol = 2)
logRegs2 <- CIN_LogReg + CLE_LogReg + DET_LogReg + JAX_LogReg + plot_layout(ncol = 2)
logRegs3 <- MTL_LogReg + NY_LogReg + NOR_LogReg + PHL_LogReg + plot_layout(ncol = 2)
logRegs4 <- PIT_LogReg + TOR_LogReg + WDC_LogReg + plot_layout(ncol = 2)

ggsave(filename = "analysis/figures/test.pdf", plot = logRegs2, device = "pdf", 
       width = 8, height = 8, units = "in", dpi = 300) 

##############################################
#### ANALYSIS OF CLINES ACROSS ALL CITIES ####
##############################################

## POPULATION-MEAN HCN FREQUENCIES ###

# Load in data with population-mean HCN frequencies
datPops <- read_csv("data-clean/AllCities_AllPopulations.csv")

# Anova for overall effect of distance, city, and whether the strenth of clines
# varies across cities (i.e., distance x city interaction)
clinesAllCities <- lm(freqHCN ~ std_distance*City, data = datPops)
summary(clinesAllCities)
Anova(clinesAllCities, type = 3)

## INDIVIDUAL PLANT PHENOTYPE DATA ##

# Read in individual plant-level data. Append std_distance
datPlants <- read_csv("data-clean/AllCities_AllPlants.csv") %>% 
  left_join(., datPops %>%select(City, Population, std_distance), by = c("City", "Population"))

# Perform logistic regression
logClinesAllCities <- glm(HCN_Result ~ std_distance*City, family = "binomial", data = datPlants)
summary(logClinesAllCities)
anova(logClinesAllCities, test = "LRT")

intercept <- coef(logClinesAllCities)["(Intercept)"]
distance_effect <- coef(logClinesAllCities)["std_distance"]

probCyan_urban <- exp(intercept + distance_effect * 0) / (1 + (exp(intercept + distance_effect * 0)))
probCyan_rural <- exp(intercept + distance_effect * 1) / (1 + (exp(intercept + distance_effect * 1)))
(probCyan_rural - probCyan_urban) / probCyan_urban # Relative increase in probability of being cyanogenic from rural --> urban


# Same model as above but specified using different syntax
logClinesAllCitiesFreqs <- glm(freqHCN ~ std_distance*City, family = "binomial", data = datPops, weights = n_HCN)
summary(logClinesAllCitiesFreqs)
anova(logClinesAllCitiesFreqs, test = "LRT")


###############################################
#### CORRELATIONS AMONG CLIMATIC VARIABLES ####
###############################################

# Load in city summary dataset
citySummaryData <- read_csv("data-clean/citySummaryData.csv") 

# Select weather variables
weather_data <- citySummaryData %>%
  select(Latitude, Longitude, annualAI, monthlyPET, annualPET,
         monthlyPrecip, mwtBio, mstBio, smd, snow_depth, snowfall,
         daysNegNoSnow) %>%
  as.matrix()

# Create correlation matrix
corr_mat_object <- rcorr(weather_data, type = "pearson")
corr_mat <- round(corr_mat_object$r, 3)

# Assign P-values to upper triangle of correlation matrix
corr_mat[upper.tri(corr_mat)] <- round(corr_mat_object$P[upper.tri(corr_mat_object$P)], 3)

# Write correlation matrix to disk
corr_mat <- as.data.frame(corr_mat) %>%
  rownames_to_column()

##################################################
#### ANALYSIS PREDICTING MEAN HCN FREQUENCIES ####
##################################################

# Get range of HCN frequencies
citySummaryData %>%
  arrange(freqHCN) %>%
  select(City, freqHCN)

#Run Models for each environmental variable predicting the strength of clines
summary(lm(freqHCN ~ Latitude, data = citySummaryData)) # KEEP
summary(lm(freqHCN ~ Longitude, data = citySummaryData)) # REMOVE
summary(lm(freqHCN ~ annualAI, data = citySummaryData)) # MARNGINAL. KEEP
summary(lm(freqHCN ~ monthlyPET, data = citySummaryData)) # KEEP
summary(lm(freqHCN ~ annualPET, data = citySummaryData)) # KEEP
summary(lm(freqHCN ~ monthlyPrecip, data = citySummaryData)) # KEEP
summary(lm(freqHCN ~ mwtBio, data = citySummaryData)) # KEEP
summary(lm(freqHCN ~ mstBio, data = citySummaryData)) # KEEP
summary(lm(freqHCN ~ smd, data = citySummaryData)) # KEEP
summary(lm(freqHCN ~ snow_depth, data = citySummaryData)) # REMOVE
summary(lm(freqHCN ~ snowfall, data = citySummaryData)) # KEEP
summary(lm(freqHCN ~ daysNegNoSnow, data = citySummaryData)) # KEEP

# Predictors kept if P <= 0.1 from above models
envPredictorsHCN <- citySummaryData %>%
  column_to_rownames('City') %>%
  select(annualAI, annualPET,
         monthlyPET, monthlyPrecip,
         mwtBio, mstBio, smd, snowfall,
         daysNegNoSnow)

# AnnualAI, daysNegNoSnow, and smd show only moderate correlations with other
# variables. These will be kept as distinct predictors. monthlyPET will be
# eliminated since it is highly correlated with, and measures the same thing
# as, annualPET. The remaining variables all show strong correlations and will
# be reduced through PCA. See text S1 for details and table SX for correlation matrix

envPredictorsHCN_forPCA <- envPredictorsHCN %>%
  select(annualPET, monthlyPrecip, mwtBio, mstBio,
         snowfall)

# Perform PCA of remaining environmental variables 
envPCAHCN <- prcomp(envPredictorsHCN_forPCA, center = TRUE, scale = TRUE)

envPCAHCN_eigen <- get_eigenvalue(envPCAHCN)
envPCAHCN_eigen # Percent variation and cummulative % of each PC

envPCAHCN_vars <- get_pca_var(envPCAHCN)
envPCAHCN_vars$coord  # Variable loadings
envPCAHCN_vars$contrib  # Variable contributions to the PCs

envPCAHCN_inds <- get_pca_ind(envPCAHCN)
envPCAHCN_inds$coord  # City loadings
envPCAHCN_inds$contrib  # City contributions to the PCs

# Assess how many PCs to keep based on broken stick method
# Keep PCs if 'Inertia' is above broken stick line
screeplot(envPCAHCN, bstick = TRUE)

# Pull loadings out for each city and merge PC1 (92.8% variation explained) with 
# analysis dataset
citySummaryData <- envPCAHCN_inds$coord  %>% # PCA scores for cities. Can be extracted and used in regression
  as.data.frame() %>%
  select(Dim.1) %>%
  rename(PC1_HCN = Dim.1) %>%
  rownames_to_column("City") %>%
  merge(., citySummaryData, by = "City")

# citySummaryData$PC1_HCN_squared <- citySummaryData$PC1_HCN^2
# Model for change in mean HCN frequency
HCNfreqMod <- lm(freqHCN ~ PC1_HCN + annualAI + daysNegNoSnow +
                   smd, 
                 data = citySummaryData)
summary(HCNfreqMod)

# Model selection and averaging of mena HCN frequency model
options(na.action = "na.fail")
HCNfreqMod_dredge <- dredge(HCNfreqMod, rank = "AICc", 
                            evaluate = TRUE,
                            extra = c("R^2", "adjR^2", F = function(x)
                              summary(x)$fstatistic[[1]]))
options(na.action = "na.omit")

# Get all dredge models as data frame and summarise output
HCN_dredge_models <- as.data.frame(HCNfreqMod_dredge)
models_HCN <- get.models(HCNfreqMod_dredge, subset = delta < 2)
HCN_modAvg <- model.avg(models_HCN)
summary(HCN_modAvg)

#######################################################
#### ANALYSIS OF FACTORS PREDICTING CLINE STRENGTH ####
#######################################################

## POPULATION-MEAN DATA ##

# Generate reduced dataset that excludes Tampa, which is fixed for HCN
citySummaryDataForAnalysis <- citySummaryData %>%
  filter(City != "Tampa")

#Run Models for each environmental variable predicting the strength of clines
summary(lm(betaLinOnly ~ Latitude, data = citySummaryDataForAnalysis)) # KEEP
summary(lm(betaLinOnly ~ Longitude, data = citySummaryDataForAnalysis)) # REMOVE
summary(lm(betaLinOnly ~ annualAI, data = citySummaryDataForAnalysis)) # REMOVE
summary(lm(betaLinOnly ~ monthlyPET, data = citySummaryDataForAnalysis)) # REMOVE
summary(lm(betaLinOnly ~ annualPET, data = citySummaryDataForAnalysis)) # REMOVE
summary(lm(betaLinOnly ~ monthlyPrecip, data = citySummaryDataForAnalysis)) # REMOVE
summary(lm(betaLinOnly ~ mwtBio, data = citySummaryDataForAnalysis)) # KEEP
summary(lm(betaLinOnly ~ mstBio, data = citySummaryDataForAnalysis)) # MARGINAL. KEEP.
summary(lm(betaLinOnly ~ smd, data = citySummaryDataForAnalysis)) # REMOVE
summary(lm(betaLinOnly ~ snow_depth, data = citySummaryDataForAnalysis)) # MARGINAL. KEEP.
summary(lm(betaLinOnly ~ snowfall, data = citySummaryDataForAnalysis)) # KEEP
summary(lm(betaLinOnly ~ daysNegNoSnow, data = citySummaryDataForAnalysis)) # REMOVE

# Pull out columns with remaining environmental predictors.
# Predictors kept if P <= 0.1 from above models
envPredictorsSlope <- citySummaryDataForAnalysis %>%
  column_to_rownames('City') %>%
  select(mwtBio, mstBio, snow_depth, snowfall)

# Perform PCA of remaining environmental variables 
envPCAslope <- prcomp(envPredictorsSlope, center = TRUE, scale = TRUE)

envPCAslope_eigen <- get_eigenvalue(envPCAslope)
envPCAslope_eigen # Percent variation and cummulative % of each PC

envPCAslope_vars <- get_pca_var(envPCAslope)
envPCAslope_vars$coord  # Variable loadings
envPCAslope_vars$contrib  # Variable contributions to the PCs

envPCAslope_inds <- get_pca_ind(envPCAslope)
envPCAslope_inds$coord  # City loadings
envPCAslope_inds$contrib  # City contributions to the PCs

# Assess how many PCs to keep based on broken stick method
# Keep PCs if 'Inertia' is above broken stick line
screeplot(envPCAslope, bstick = TRUE)

# Pull loadings out for each city and merge PC1 (92.8% variation explained) with 
# analysis dataset
citySummaryDataForAnalysis <- envPCAslope_inds$coord  %>% # PCA scores for cities. Can be extracted and used in regression
  as.data.frame() %>%
  select(Dim.1) %>%
  rename(PC1_Slope = Dim.1) %>%
  rownames_to_column("City") %>%
  merge(., citySummaryDataForAnalysis, by = "City")

SlopeMod <- lm(betaLinOnly ~ PC1_Slope, data = citySummaryDataForAnalysis)
summary(SlopeMod)

## LOGISTIC REGRESSION DATA ##

#Run Models for each environmental variable predicting the strength of clines
summary(lm(betaLog ~ Latitude, data = citySummaryDataForAnalysis)) # KEEP
summary(lm(betaLog ~ Longitude, data = citySummaryDataForAnalysis)) # REMOVE
summary(lm(betaLog ~ annualAI, data = citySummaryDataForAnalysis)) # REMOVE
summary(lm(betaLog ~ monthlyPET, data = citySummaryDataForAnalysis)) # REMOVE
summary(lm(betaLog ~ annualPET, data = citySummaryDataForAnalysis)) # KEEP
summary(lm(betaLog ~ monthlyPrecip, data = citySummaryDataForAnalysis)) # KEEP
summary(lm(betaLog ~ mwtBio, data = citySummaryDataForAnalysis)) # KEEP
summary(lm(betaLog ~ mstBio, data = citySummaryDataForAnalysis)) # MARGINAL. KEEP.
summary(lm(betaLog ~ smd, data = citySummaryDataForAnalysis)) # KEEP
summary(lm(betaLog ~ snow_depth, data = citySummaryDataForAnalysis)) # MARGINAL. KEEP.
summary(lm(betaLog ~ snowfall, data = citySummaryDataForAnalysis)) # KEEP
summary(lm(betaLog ~ daysNegNoSnow, data = citySummaryDataForAnalysis)) # REMOVE

# Pull out columns with remaining environmental predictors.
# Predictors kept if P <= 0.1 from above models
envPredictorsSlopeLog <- citySummaryDataForAnalysis %>%
  column_to_rownames('City') %>%
  select(annualPET, monthlyPrecip, mwtBio, mstBio, snow_depth, snowfall)

# Perform PCA of remaining environmental variables 
envPCAslopeLog <- prcomp(envPredictorsSlopeLog, center = TRUE, scale = TRUE)

envPCAslope_eigenLog <- get_eigenvalue(envPCAslopeLog)
envPCAslope_eigenLog # Percent variation and cummulative % of each PC

envPCAslope_varsLog <- get_pca_var(envPCAslopeLog)
envPCAslope_varsLog$coord  # Variable loadings
envPCAslope_varsLog$contrib  # Variable contributions to the PCs

envPCAslope_indsLog <- get_pca_ind(envPCAslopeLog)
envPCAslope_indsLog$coord  # City loadings
envPCAslope_indsLog$contrib  # City contributions to the PCs

# Assess how many PCs to keep based on broken stick method
# Keep PCs if 'Inertia' is above broken stick line
screeplot(envPCAslopeLog, bstick = TRUE)

# Pull loadings out for each city and merge PC1 (92.8% variation explained) with 
# analysis dataset
citySummaryDataForAnalysis <- envPCAslope_indsLog$coord  %>% # PCA scores for cities. Can be extracted and used in regression
  as.data.frame() %>%
  select(Dim.1) %>%
  rename(PC1_Slope = Dim.1) %>%
  rownames_to_column("City") %>%
  merge(., citySummaryDataForAnalysis, by = "City")

SlopeModLog <- lm(betaLog ~ PC1_Slope, data = citySummaryDataForAnalysis)
summary(SlopeModLog)

#### ANALYSIS OF HAPLOTYPES BY HABITAT TYPE ####

# Load in data
haplotype_data <- read_csv("data-clean/haplotypeData.csv")

# Number of plants per population and city
haplotype_data %>% 
  group_by(City, Habitat, haplotype_Ac) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  group_by(City, Habitat) %>% 
  summarise(sum = sum(count))

## Ac LOCUS ##

# Counts and proportion of different haplotypes at Ac Locus
haplo_counts_Ac <- haplotype_data %>%
  group_by(haplotype_Ac) %>%
  summarise(count = n(),
            prop = round(count / nrow(haplotype_data), 3))

haplo_counts_Ac_City <- haplotype_data %>%
  group_by(haplotype_Ac, City) %>%
  summarise(count = n(),
            prop = round(count / nrow(haplotype_data), 3))

# Correct/incorrect calls at Ac locus
haplo_validation_Ac <- haplo_counts_Ac %>%
  mutate(validation = case_when(
    grepl("Ac", haplotype_Ac) ~ "Glucoside_positive",
    grepl("unk", haplotype_Ac) ~ "Invalid",
    TRUE ~ "Valid"
  )) %>%
  group_by(validation) %>%
  summarise(total_prop = sum(prop))

# Number of "ac" haplotypes by city and habitat
richnessHaplotAc <- haplotype_data %>% 
  filter(haplotype_Ac != "unk" & haplotype_Ac != "Ac") %>% 
  group_by(City, Habitat) %>% 
  distinct_at(., "haplotype_Ac") %>% 
  summarise(richness = n())

# Model testing for variation in "ac" haplotype richness by habitat
AcLocusMod_Rich <- aov(richness ~ Habitat, data = richnessHaplotAc)
summary(AcLocusMod_Rich)

## Li LOCUS ##

# Counts and proportion of different haplotypes at Li Locus
haplo_counts_Li <- haplotype_data %>%
  group_by(haplotype_Li) %>%
  summarise(count = n(),
            prop = round(count / nrow(haplotype_data), 3))

haplo_counts_Li_City <- haplotype_data %>%
  group_by(haplotype_Li, City) %>%
  summarise(count = n(),
            prop = round(count / nrow(haplotype_data), 3))

# Correct/incorrect calls at Li locus
haplo_validation_Li <- haplo_counts_Li %>%
  mutate(validation = case_when(
    grepl("Li", haplotype_Li) ~ "Enzyme_positive",
    grepl("unk", haplotype_Li) ~ "Invalid",
    TRUE ~ "Valid"
  )) %>%
  group_by(validation) %>%
  summarise(total_prop = sum(prop))

# Number of "li" haplotypes by city and habitat
richnessHaplotLi <- haplotype_data %>% 
  filter(haplotype_Li != "unk" & haplotype_Li != "Li") %>% 
  group_by(City, Habitat) %>% 
  distinct_at(., "haplotype_Li") %>% 
  summarise(richness = n())
  
# Model testing for variation in "ac" haplotype richness by habitat
LiLocusMod_Rich <- aov(richness ~ Habitat, data = richnessHaplotLi)
summary(LiLocusMod_Rich)

#### TABLES ####

## TABLE 1 

# Linear cline models only. Part of table 1
linearClines <- linearClineModelOnly(city_dataframes)

# Extract number of populations and plants per city. Merge with linear clines model output above
table1 <- read.csv("data-clean/AllCities_AllPlants.csv") %>%
  group_by(City) %>%
  summarise(numPops = n_distinct(Population),
            numPlants = n()) %>%
  merge(., linearClines, by = "City")

# Write table 1 to disk
write_csv(table1, "analysis/tables/main-text/Table1_cityClineSummary.csv")

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

write_csv(corr_mat, path = "analysis/tables/supplemental/TableS2_weatherCorrMat.csv", 
          col_names = TRUE)

## TABLE S3

# Written to disk by `writeClineResults` function

## TABLE S4

write_csv(HCN_dredge_models, path = "analysis/tables/supplemental/TableS4_HCN_dredge_output.csv", 
          col_names = TRUE)

## TABLE S5

tableS5 <- as.data.frame(rbind(
  c("Full"),
  summary(HCN_modAvg)$coefmat.full)) %>%
  rownames_to_column()

write_csv(tableS5, "analysis/tables/supplemental/TableS5_freqHCN_FullModelAvg.csv")

#### FIGURES ####

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
  do(mod = lm(freqHCN ~ std_distance, data = .)) %>%
  tidy(., mod) %>%
  mutate(intercept = estimate[term == "(Intercept)"]) %>%
  filter(term == "std_distance") %>%
  mutate(predicted = intercept + (estimate * 1.0)) %>%
  select(City, p.value, predicted) %>%
  merge(., datPops, by = "City", all.y = TRUE) %>%
  mutate(significant = ifelse(p.value < 0.05, "Yes", "No")) %>%
  ggplot(., aes(x = std_distance, y = freqHCN)) +  
    geom_line(stat = "smooth", method="lm", aes(linetype = significant,
                                                group = City), 
              alpha = 0.7, size = 1) +
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

ggsave(filename = "analysis/figures/main-text/from_R/Figure2_HCN-by-distance.pdf", 
       plot = HCN_by_city, device = 'pdf', units = 'in',
       width = 12, height = 8, dpi = 600)

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

ggsave(filename = "analysis/figures/main-text/from_R/Figure3a_HCN-by-NumDaysNegNoSnow.pdf", 
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

ggsave(filename = "analysis/figures/main-text/from_R/Figure3b_HCN-by-PC1.pdf", 
       plot = HCN_by_PC1, device = 'pdf', units = 'in',
       width = 5, height = 5, dpi = 600)


## FIGURE 4 

citySummaryDataForAnalysis <- citySummaryDataForAnalysis %>% 
  left_join(., citySummaryData %>% select(City, abbr), by = "City")

# Figure 4. Slope against PC1
Slope_by_PC1 <- citySummaryDataForAnalysis %>%
  ggplot(., aes(x = PC1_Slope, y = cyanSlopeForAnalysis)) +
  # geom_point(size = 2.5) +
  geom_smooth(method = "lm", size = 1.5, colour = "black", 
              se = FALSE) +
  # scale_x_continuous(breaks = seq(from = 0, to = 35, by = 5)) +
  xlab("PC1 (92.8%)") + ylab("Slope of HCN cline") +
  geom_text(aes(label = abbr), vjust = 0, hjust = 0) + 
  ng1
Slope_by_PC1

ggsave(filename = "analysis/figures/main-text/from_R/Figure4_Slope-by-PC1.pdf", 
       plot = Slope_by_PC1, device = 'pdf', units = 'in',
       width = 5, height = 5, dpi = 600)

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

ggsave(filename = "analysis/figures/supplemental/from_R/FigureS2_HCN_by_Lat.pdf", 
       plot = plotHCN_by_lat, device = "pdf", 
       width = 5, height = 5, dpi = 300)

## FIGURE S3

Ac_Haplo_Freqs <- ggplot(freqHaploAc, aes(x = Habitat, y = freqHaploAc, fill = haplotype_Ac)) +
  geom_bar(stat='identity') + 
  xlab("Habitat") + ylab("Haplotype frequency") + 
  facet_grid(~City) +
  scale_fill_manual(values = c("#F98400", "#5BBCD6")) +
  ng1 + theme(legend.position = "top",
              legend.direction = "horizontal",
              aspect.ratio=3.0, legend.text=element_text(size=10),
              legend.title = element_text(size = 0), 
              legend.key.size = unit(0.5, "cm"),
              axis.text.x = element_text(angle = 45, hjust = 1))
Ac_Haplo_Freqs

ggsave(filename = 'analysis/figures/supplemental/from_R/FigureS3_Ac_haplotype_freqs.pdf', 
       plot = Ac_Haplo_Freqs, device = "pdf", 
       width = 8, height = 5, dpi = 300)

## FIGURE S4

Li_Haplo_Freqs <- ggplot(freqHaploLi, aes(x = Habitat, y = freqHaploLi, fill = haplotype_Li)) +
  geom_bar(stat='identity') + 
  xlab("Habitat") + ylab("Haplotype frequency") + 
  facet_grid(~City) +
  scale_fill_manual(values = c("#00A08A", "#F2AD00", "#F98400", "#5BBCD6")) +
  ng1 + theme(legend.position = "top",
              legend.direction = "horizontal",
              aspect.ratio=3.0, legend.text=element_text(size=10),
              legend.title = element_text(size = 0), 
              legend.key.size = unit(0.5, "cm"),
              axis.text.x = element_text(angle = 45, hjust = 1))
Li_Haplo_Freqs

ggsave(filename = 'analysis/figures/supplemental/from_R/FigureS4_Li_haplotype_freqs.pdf', 
       plot = Li_Haplo_Freqs, device = "pdf", 
       width = 8, height = 5, dpi = 300)

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
