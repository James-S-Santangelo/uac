

### SETUP ###

# Clear environement, if necessary
rm(list=ls())
.rs.restartR()

# Load required packages
library(tidyverse)
library(broom)
library(FactoMineR)
library(factoextra)
library(MuMIn)
# library(RColorBrewer)

# Load data with population-level data for all cities
datPops <- read.csv("data-clean/AllCities_AllPopulations.csv")

#### FUNCTIONS ####

getBestFitClineModel <- function(response, data) {
  #' Assesses whether a linear or quadratic model is best fit
  #'
  #' Quadratic model is assumed to be better fit if it improves model AIC
  #' by more than 2 points.
  #' 
  #' @param response Response variable to use in cline models. One of "freqHCN",
  #' "AcHWE", or "LiHWE".
  #' @param data Dataset to use for running cline model.
  #' @return list with first element being the lm model object for the best fit model and
  #' the second being a string with the model order ("linear" or "quadratic").
  
  # Pull variables for linear model from dataset
  response <- data[, response]
  std_distance <- data[, "std_distance"]
  std_distance_squared <- data[, "std_distance_squared"]
  
  # Define linear and quadratic models for analysis of HCN with distance
  quadratic_model = lm(response ~ std_distance + std_distance_squared) # Specify quadratic model
  linear_model = update(quadratic_model, ~ . - std_distance_squared) # Specify linear model
  
  AIC_quad = AIC(quadratic_model) # Get AIC of quadratic model
  AIC_lin = AIC(linear_model) # Get AIC of linear model
  
  if (abs(AIC_quad) - abs(AIC_lin) > 2) {
    # If quadratic model AIC is > 2 from linear model AIC
    getBestFitClineModel = quadratic_model # Then best fit model is quadratic
    order = "quadratic"
  } else {
    # Otherwise (i.e. quadratic model is not better fit)
    getBestFitClineModel = linear_model # Then best fit model is linear
    order = "linear"
  }
  
  # Return list with best fit model and model order as string
  return(list(model_output = getBestFitClineModel,
              model_order = order))
}

getClineModelOutputFromOrder <-
  function(best_model_output,
           best_model_order,
           clines_results) {
    #' Retrieves the slopes and P-values for linear and quadratic (if present) terms from lm object
    #'
    #' @param best_model_output lm object representing best fit cline model
    #' @param best_model_order String representing order of best fit model
    #' @param clines_results Vector storing coefficients from model output
    #' @return clines_results: Vector storing coefficients from model output
    
    
    # Tidy model output
    tidy_output <- broom::tidy(best_model_output)
    
    # Get model output based on model order
    if (best_model_order == "quadratic") {
      # print("quadratic")
      coef_lin = round(tidy_output[tidy_output$term == "std_distance", "estimate"], 4)
      # print(coef_lin)
      pval_lin = round(tidy_output[tidy_output$term == "std_distance", "p.value"], 4)
      coef_quad = round(tidy_output[tidy_output$term == "std_distance_squared", "estimate"], 4)
      pval_quad = round(tidy_output[tidy_output$term == "std_distance_squared", "p.value"], 4)
      results <- c(coef_lin,
                   pval_lin,
                   coef_quad,
                   pval_quad)
      clines_results <-
        append(clines_results, results, after = length(clines_results))
    } else if (best_model_order == "linear") {
      # print("linear")
      coef_lin = round(tidy_output[tidy_output$term == "std_distance", "estimate"], 4)
      # print(coef_lin)
      pval_lin = round(tidy_output[tidy_output$term == "std_distance", "p.value"], 4)
      coef_quad = "NA" # NA since no quadratic term in model
      pval_quad = "NA" # NA since no quadratic term in model
      results <- c(coef_lin,
                   pval_lin,
                   coef_quad,
                   pval_quad)
      clines_results <-
        append(clines_results, results, after = length(clines_results))
    }
    return(clines_results)
  }

writeClineResults <- function(dataframe_list) {
  #' Writes cline model output for HCN, Ac, and Li to table
  #'
  #' @param dataframe_list A list containing individual dataframes with
  #' the frequency of HCN, Ac (if present), and Li (if present), and
  #' standardized distance to the urban centre and distance squared
  #' (for quadratic cline model)
  #' @return None: Writes dataframe to disk and assigns to global environment
  
  # Initialize dataset that will hold model outputs
  modelOutputData <- data.frame(
    City = character(),
    cyanSlopeLin = numeric(),
    pvalCyanSlopeLin = numeric(),
    cyanSlopeQuad = numeric(),
    pvalSlopeQuad = numeric(),
    AcSlopeLin = numeric(),
    pvalAcSlopeLin = numeric(),
    AcSlopeQuad = numeric(),
    pvalAcSlopeQuad = numeric(),
    LiSlopeLin = numeric(),
    pvalLoSlopeLin = numeric(),
    LiSlopeQuad = numeric(),
    pvalLiSlopeQuad = numeric(),
    stringsAsFactors = FALSE
  )
  
  clineOrderData <- data.frame(
    City = character(),
    cyanOrder = character(),
    AcOrder = character(),
    LiOrder = character(),
    stringsAsFactors = FALSE
  )
  
  # Iterate over list containing dataframes
  for (i in 1:length(dataframe_list)) {
    # Retrieve dataframe from list
    dataframe = dataframe_list[[i]]
    
    # Extract city as character
    city = as.character(unique(dataframe$City))
    
    # Initialize empty vectors that will be grown and appendeded to written dataframes
    clines_results <- c()
    order_results <- c()
    
    ## HCN ##
    
    # Identify best fit model, return model and string with model order
    best_model_HCN <-
      getBestFitClineModel("freqHCN", dataframe)
    best_model_output_HCN <- best_model_HCN[["model_output"]]
    best_model_order_HCN <- best_model_HCN[["model_order"]]
    
    # Append output for HCN cline to cline results vector.
    clines_results <-
      getClineModelOutputFromOrder(best_model_output_HCN,
                                   best_model_order_HCN,
                                   clines_results)
    # print(clines_results)
    
    # Append model order to order results
    order_results <-
      append(order_results,
             best_model_order_HCN,
             after = length(order_results))
    
    # Write best fit HCN model to disk
    write_csv(
      tidy(best_model_output_HCN),
      path = sprintf(
        "analysis/inividual-cline-models/HCN/%s-HCN-cline-model.csv",
        city
      ),
      append = FALSE
    )
    
    ## Ac and Li, IF PRESENT ##
    
    # If any row contains values for inferred Ac and Li HWE allele frequencies
    if (any(!is.na(dataframe[dataframe$City == city, c("AcHWE", "LiHWE")])) == FALSE) {
      clines_results <-
        append(
          clines_results,
          c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"),
          after = length(clines_results)
        )
      order_results <-
        append(order_results, c("NA", "NA"), after = length(order_results))
    } else {
      ## AC ##
      
      # Identify best fit model, return model and string with model order
      best_model_Ac <-
        getBestFitClineModel("AcHWE", dataframe)
      best_model_output_Ac <- best_model_Ac[["model_output"]]
      best_model_order_Ac <- best_model_Ac[["model_order"]]
      
      clines_results <-
        getClineModelOutputFromOrder(best_model_output_Ac,
                                     best_model_order_Ac,
                                     clines_results)
      order_results <-
        append(order_results,
               best_model_order_Ac,
               after = length(order_results))
      
      # Write best fit Ac model to disk
      write_csv(
        tidy(best_model_output_Ac),
        path = sprintf(
          "analysis/inividual-cline-models/Ac/%s-Ac-cline-model.csv",
          city
        ),
        append = FALSE
      )
      
      ## LI ##
      
      # Identify best fit model, return model and string with model order
      best_model_Li <-
        getBestFitClineModel("LiHWE", dataframe)
      best_model_output_Li <- best_model_Li[["model_output"]]
      best_model_order_Li <- best_model_Li[["model_order"]]
      
      clines_results <-
        getClineModelOutputFromOrder(best_model_output_Li,
                                     best_model_order_Li,
                                     clines_results)
      # print(clines_results)
      order_results <-
        append(order_results,
               best_model_order_Li,
               after = length(order_results))
      
      # Write best fit Li model to disk
      write_csv(
        tidy(best_model_output_Li),
        path = sprintf(
          "analysis/inividual-cline-models/Li/%s-Li-cline-model.csv",
          city
        ),
        append = FALSE
      )
    }
    
    # Append results vector to dataframe
    modelOutputData[i, ] <- c(city, clines_results)
    clineOrderData[i, ] <-  c(city, order_results)
  }
  # Assign dataframe to global environment and write to disk
  assign("clineModelOutput", modelOutputData, envir = .GlobalEnv)
  write_csv(
    modelOutputData,
    "analysis/clineModelOutput.csv",
    na = "NA",
    append = FALSE
  )
  
  assign("clineModelOrder", clineOrderData, envir = .GlobalEnv)
  write_csv(clineOrderData,
            "analysis/clineOrderData.csv",
            na = "NA",
            append = FALSE)
}

##Function for changing upper panel in 'pairs' correlation matrix to show correlation
#coefficients and p-values
panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y, use = "complete.obs", method = "pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.001) txt2 <- paste("p= ", "<0.001", sep = "")
  text(0.5, 0.4, txt2)
}

#Function to add least squares regression line to lower panel of 'pairs' scatterplot matrix
lsline = function(x,y) {
  points(x,y,pch=".")
  abline(lsfit(x,y),col="blue")
}

#### OUPUT FROM INDIVIDUAL CLINES ####

# Split cities into different dataframes stored as list
city_dataframes <- split(datPops, datPops$City)

# Write cline results to disk and global environment
writeClineResults(city_dataframes)


#### ANALYSIS OF CLINES ACROSS ALL CITIES ####

## HCN ##

clinesAllCities <- lm(freqHCN ~ std_distance*City, data = datPops)
summary(clinesAllCities)
car::Anova(clinesAllCities, type = 3)

# colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', 
#             '#f58231', '#911eb4', '#46f0f0', '#f032e6', 
#             '#bcf60c', '#fabebe', '#008080', '#e6beff', 
#             '#9a6324', '#000075', '#800000', '#aaffc3')

## AC ##

datIndAlleles_Ac <- datPops %>%
  filter(!is.na(AcHWE))
clinesAllCities_Ac <- lm(AcHWE ~ std_distance*City, data = datIndAlleles_Ac)
summary(clinesAllCities_Ac)
car::Anova(clinesAllCities_Ac, type = 3)

## LI ##

datIndAlleles_Li <- datPops %>%
  filter(!is.na(LiHWE))
clinesAllCities_Ac <- lm(LiHWE ~ std_distance*City, data = datIndAlleles_Li)
summary(datIndAlleles_Li)
car::Anova(clinesAllCities_Ac, type = 3)

#### ANALYSIS OF FACTORS PREDICTING CLINE STRENGTH ####

# Load in city summary dataset
citySummaryData <- read_csv("data-clean/citySummaryData.csv")

# Generate reduced dataset that excludes Tampa, which is fixed for HCN
citySummaryDataForAnalysis <- citySummaryData %>%
  filter(City != "Tampa")

## CORRELATIONS AMOMG PREDICTORS

# Pairwise correlation matrix
pairs(envPredictors, upper.panel = panel.cor, lower.panel = lsline)

#Run Models for each environmental variable predicting the strength of clines
summary(lm(cyanSlopeForAnalysis ~ Latitude, data = citySummaryDataForAnalysis)) # KEEP
summary(lm(cyanSlopeForAnalysis ~ Longitude, data = citySummaryDataForAnalysis)) # REMOVE
summary(lm(cyanSlopeForAnalysis ~ annualAI, data = citySummaryDataForAnalysis)) # REMOVE
summary(lm(cyanSlopeForAnalysis ~ monthlyPET, data = citySummaryDataForAnalysis)) # REMOVE
summary(lm(cyanSlopeForAnalysis ~ annualPET, data = citySummaryDataForAnalysis)) # REMOVE
summary(lm(cyanSlopeForAnalysis ~ monthlyPrecip, data = citySummaryDataForAnalysis)) # REMOVE
summary(lm(cyanSlopeForAnalysis ~ mwtBio, data = citySummaryDataForAnalysis)) # KEEP
summary(lm(cyanSlopeForAnalysis ~ mstBio, data = citySummaryDataForAnalysis)) # MARGINAL. KEEP.
summary(lm(cyanSlopeForAnalysis ~ smd, data = citySummaryDataForAnalysis)) # REMOVE
summary(lm(cyanSlopeForAnalysis ~ snow_depth, data = citySummaryDataForAnalysis)) # MARGINAL. KEEP.
summary(lm(cyanSlopeForAnalysis ~ snowfall, data = citySummaryDataForAnalysis)) # KEEP
summary(lm(cyanSlopeForAnalysis ~ mwtWea, data = citySummaryDataForAnalysis)) # REMOVE
summary(lm(cyanSlopeForAnalysis ~ daysNegNoSnow, data = citySummaryDataForAnalysis)) # REMOVE

# Pull out columns with remaining environmental predictors.
# Predictors kept if P <= 0.1 from above models
envPredictorsSlope <- citySummaryDataForAnalysis %>%
  column_to_rownames('City') %>%
  select(Latitude, mwtBio, mstBio, 
         snow_depth, snowfall)

# Visualize correlations among remaining predictors.
pairs(envPredictorsSlope, upper.panel = panel.cor, lower.panel = lsline)

# Perform PCA of remaining environmental variables due to high correlations
envPCAslope <- PCA(envPredictorsSlope, scale.unit = T, graph = F)
envPCAslope$eig # Eigenvalues, percent variation and cummulative percent variation
envPCAslope$var$coord # Variable loadings

# Pull loadings out for each city and merge PC1 (92.8% variation explained) with 
# analysis dataset
citySummaryDataForAnalysis <- envPCAslope$ind$coord  %>% # PCA scores for cities. Can be extracted and used in regression
  as.data.frame() %>%
  select(Dim.1) %>%
  rename(PC1 = Dim.1) %>%
  rownames_to_column("City") %>%
  merge(., citySummaryDataForAnalysis, by = "City")

# Run model with PC1 predicting the strength of clines.
summary(lm(cyanSlopeForAnalysis ~ PC1, data = citySummaryDataForAnalysis)) 




fviz_pca_biplot(envPCA,
                repel = T,
                title = "",
                col.var = "black") +
  # scale_color_manual(values = c("black", "red")) +
  ylab("PC2") + xlab("PC1") +
  ng1

citySummaryDataForAnalysis %>%
  ggplot(., aes(x = PC1, y = cyanSlopeForAnalysis)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ng1

datPops %>%
  filter(City == "Charlotte") %>%
  ggplot(., aes(x = std_distance, y = freqHCN)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2) ) +
  ng1

#### FIGURES ####

#Theme used to plot figures throughout script. Modify depending on usage (e.g. poster, talk, etc.)
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

# Plot of HCN frequency against distance for each city
# Solid line if significan. Thick black line is regression across all cities.
datPops %>%
  group_by(City) %>%
  do(mod = lm(freqHCN ~ std_distance, data = .)) %>%
  tidy(., mod) %>%
  filter(term == "std_distance") %>%
  select(City, p.value) %>%
  merge(., datPops, by = "City", all.y = TRUE) %>%
  mutate(significant = ifelse(p.value < 0.05, "Yes", "No")) %>%
  ggplot(., aes(x = std_distance, y = freqHCN)) +  
    geom_line(stat = "smooth", method="lm", aes(group = City, linetype = significant), alpha = 0.5, size = 1) +
    geom_line(stat = "smooth", method="lm", colour = "black", size = 2.5) +
    # scale_colour_manual(values = colors) +
    xlab("Standardized distance") + ylab("Frequency of HCN") + 
    scale_linetype_manual(values=c("dashed", "solid")) +
    scale_y_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
    coord_cartesian(ylim = c(0.1, 1.025)) +
    ng1

# Plot of Ac frequency against distance for each city
# Solid line if significan. Thick black line is regression across all cities.
datIndAlleles_Ac %>%
  group_by(City) %>%
  do(mod = lm(AcHWE ~ std_distance, data = .)) %>%
  tidy(., mod) %>%
  filter(term == "std_distance") %>%
  select(City, p.value) %>%
  merge(., datIndAlleles_Ac, by = "City", all.y = TRUE) %>%
  mutate(significant = ifelse(p.value < 0.05, "Yes", "No")) %>%
  ggplot(., aes(x = std_distance, y = AcHWE)) +  
  geom_line(stat = "smooth", method="lm", aes(group = City, linetype = significant), alpha = 0.5, size = 1) +
  geom_line(stat = "smooth", method="lm", colour = "black", size = 2.5) +
  # scale_colour_manual(values = colors) +
  xlab("Standardized distance") + ylab("Frequency of Ac") + 
  scale_linetype_manual(values=c("dashed", "solid")) +
  scale_y_continuous(breaks = seq(from = 0.2, to = 1, by = 0.1)) +
  coord_cartesian(ylim = c(0.2, 1.025)) +
  ng1

# Plot of Li frequency against distance for each city
# Solid line if significan. Thick black line is regression across all cities.
datIndAlleles_Li %>%
  group_by(City) %>%
  do(mod = lm(LiHWE ~ std_distance, data = .)) %>%
  tidy(., mod) %>%
  filter(term == "std_distance") %>%
  select(City, p.value) %>%
  merge(., datIndAlleles_Li, by = "City", all.y = TRUE) %>%
  mutate(significant = ifelse(p.value < 0.05, "Yes", "No")) %>%
  ggplot(., aes(x = std_distance, y = LiHWE)) +  
  geom_line(stat = "smooth", method="lm", aes(group = City, linetype = significant), alpha = 0.5, size = 1) +
  geom_line(stat = "smooth", method="lm", colour = "black", size = 2.5) +
  # scale_colour_manual(values = colors) +
  xlab("Standardized distance") + ylab("Frequency of Li") + 
  scale_linetype_manual(values=c("dashed", "solid")) +
  scale_y_continuous(breaks = seq(from = 0.1, to = 1, by = 0.1)) +
  coord_cartesian(ylim = c(0.1, 1.025)) +
  ng1












