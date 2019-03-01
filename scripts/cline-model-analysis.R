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
library(car)
library(FactoMineR)
library(factoextra)
library(vegan)
library(RColorBrewer)
library(MuMIn)
library(Hmisc)

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
    freqHCN = character(),
    AcHWE = character(),
    LiHWE = character(),
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
    
    # If no values for inferred Ac and Li HWE allele frequencies
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

linearClineModelOnly <- function(dataframe_list){
  
  
  # Initialize dataset that will hold model outputs
  modelOutputData <- data.frame(
    City = character(),
    betaHCN = numeric(),
    pvalHCN = numeric(),
    betaAc = numeric(),
    pvalAc = numeric(),
    betaLi = numeric(),
    pvalLi = numeric(),
    stringsAsFactors = FALSE
  )
  
  
  for (i in 1:length(dataframe_list)) {
    
    # Retrieve dataframe from list
    dataframe = dataframe_list[[i]]
    
    # Extract city as character
    city = as.character(unique(dataframe$City))
    
    # Initialize vector to hold results
    cline_results <- c()
    
    # Run linear model
    clineModel <- lm(freqHCN ~ std_distance, data = dataframe)
    
    #Extract relavent coeficient
    betaHCN <- round(summary(clineModel)$coefficients["std_distance", "Estimate"], 4)
    pvalHCN <- round(summary(clineModel)$coefficients["std_distance", "Pr(>|t|)"], 4)
    
    # print(city)
    # print(betaHCN)
    # print(pvalHCN)
    
    cline_results <- append(cline_results,
                            c(betaHCN, pvalHCN),
                            after = length(cline_results))
    
    ## Ac and Li, IF PRESENT ##
    
    # If no values for inferred Ac and Li HWE allele frequencies
    if (any(!is.na(dataframe[dataframe$City == city, c("AcHWE", "LiHWE")])) == FALSE) {
      cline_results <- append(cline_results,
          c("NA", "NA", "NA", "NA"),
          after = length(cline_results)
        )
    }else{
      
      # AC
      
      clineModelAc <- lm(AcHWE ~ std_distance, data = dataframe)
      betaAc <- round(summary(clineModelAc)$coefficients["std_distance", "Estimate"], 4)
      pvalAc <- round(summary(clineModelAc)$coefficients["std_distance", "Pr(>|t|)"], 4)
      
      clineModelLi <- lm(LiHWE ~ std_distance, data = dataframe)
      betaLi <- round(summary(clineModelLi)$coefficients["std_distance", "Estimate"], 4)
      pvalLi <- round(summary(clineModelLi)$coefficients["std_distance", "Pr(>|t|)"], 4)
      
      cline_results <- append(cline_results,
                              c(betaAc, pvalAc, betaLi, pvalLi),
                              after = length(cline_results))
    }
    # print(city)
    # print(cline_results)
    modelOutputData[i, ] <- c(city, cline_results)
  }  
  return(modelOutputData)
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

#### ANALYSIS PREDICTING MEAN HCN FREQUENCIES ####

# Load in city summary dataset
citySummaryData <- read_csv("data-clean/citySummaryData.csv")

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

# Pull out columns with remaining environmental predictors.
# Predictors kept if P <= 0.1 from above models
envPredictorsHCN <- citySummaryData %>%
  column_to_rownames('City') %>%
  select(annualAI, annualPET,
         monthlyPET, monthlyPrecip,
         mwtBio, mstBio, smd, snowfall,
         daysNegNoSnow)


# AnnualAI, daysNegNoSnow, and smd show only moderate correlations with other
# variables. These will be kept as disting predictors. monthlyPET will be
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

# Get all dredge models as data frame and save to disk
HCN_dredge_models <- as.data.frame(HCNfreqMod_dredge)
write_csv(HCN_dredge_models, path = "analysis/HCN_dredge_output.csv", col_names = TRUE)

models_HCN <- get.models(HCNfreqMod_dredge, subset = delta < 2)
HCN_modAvg <- model.avg(models_HCN)
summary(HCN_modAvg)
HCN_modAvg$coefArray

#### ANALYSIS OF FACTORS PREDICTING CLINE STRENGTH ####

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
summary(lm(cyanSlopeForAnalysis ~ daysNegNoSnow, data = citySummaryDataForAnalysis)) # REMOVE

# Pull out columns with remaining environmental predictors.
# Predictors kept if P <= 0.1 from above models
envPredictorsSlope <- citySummaryDataForAnalysis %>%
  column_to_rownames('City') %>%
  select(mwtBio, mstBio, snow_depth, snowfall)

# Visualize correlations among remaining predictors.
pairs(envPredictorsSlope, upper.panel = panel.cor, lower.panel = lsline)

# All remaining variables highly correlated. Perform PCA.

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

SlopeMod <- lm(cyanSlopeForAnalysis ~ PC1_Slope, data = citySummaryDataForAnalysis)
summary(SlopeMod)

#### ANALYSIS OF HAPLOTYPES BY HABITAT TYPE ####

# Load in data
haplotype_data <- read_csv("data-clean/haplotypeData.csv")

## AC LOCUS ##

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

# Relative frequency of haplotypes at Ac Locus
freqHaploAc <- haplotype_data %>%
  group_by(City, Habitat, haplotype_Ac) %>%
  filter(haplotype_Ac != "unk" & haplotype_Ac != "Ac") %>%
  summarise(nAc = n()) %>%
  mutate(freqHaploAc = nAc / sum(nAc)) %>%
  ungroup() %>%
  mutate(haplotype_Ac = as.factor(haplotype_Ac),
         haplotype_Ac = fct_reorder(haplotype_Ac, freqHaploAc)) %>%
  complete(City, Habitat, haplotype_Ac, fill = list(freqHaploAc = 0))

# Dataframe with simpson's diversity for haplotypes
simpsDivAc <- haplotype_data %>%
  group_by(City, Habitat, haplotype_Ac) %>%
  filter(haplotype_Ac != "unk" & haplotype_Ac != "Ac") %>%
  summarize(n = n(),
            nMin1 = n - 1, 
            n_nMin1 = n * nMin1) %>%
  ungroup() %>%
  group_by(City, Habitat) %>%
  summarize(N = sum(n),
            sum_n_nMin1 = sum(n_nMin1),
            simpson = 1 - (sum_n_nMin1 / (N* (N - 1))))

# Model testing for variation in Ac haplotype relative frequency by haplotype and habitat type
AcLocusMod <- lm(freqHaploAc ~ Habitat*haplotype_Ac, data = freqHaploAc)
car::Anova(AcLocusMod, type = 3)

# Model testing for variation in Ac haplotype simpson's diversity by haplotype and habitat type
AcLocusMod_Simp <- lm(simpson ~ Habitat, data = simpsDivAc)
summary(AcLocusMod_Simp)

# Get Simpson's index in each habitat
simpsDivAc %>%
  group_by(Habitat) %>%
  summarize(meanSimp = mean(simpson))

## LI LOCUS ##

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

# Relative frequency of haplotypes at Li Locus
freqHaploLi <- haplotype_data %>%
  group_by(City, Habitat, haplotype_Li) %>%
  filter(haplotype_Li != "unk" & haplotype_Li != "Li") %>%
  summarise(nLi = n()) %>%
  mutate(freqHaploLi = nLi / sum(nLi)) %>%
  ungroup() %>%
  mutate(haplotype_Li = as.factor(haplotype_Li),
         haplotype_Li = fct_reorder(haplotype_Li, freqHaploLi)) %>%
  complete(City, Habitat, haplotype_Li, fill = list(freqHaploLi = 0))

# Dataframe with simpson's diversity for haplotypes
simpsDivLi <- haplotype_data %>%
  group_by(City, Habitat, haplotype_Li) %>%
  filter(haplotype_Li != "unk" & haplotype_Li != "Li") %>%
  summarize(n = n(),
            nMin1 = n - 1, 
            n_nMin1 = n * nMin1) %>%
  ungroup() %>%
  group_by(City, Habitat) %>%
  summarize(N = sum(n),
            sum_n_nMin1 = sum(n_nMin1),
            simpson = 1 - (sum_n_nMin1 / (N* (N - 1))))

# Model testing for variation in Li haplotype relative frequency by haplotype and habitat type
LiLocusMod <- lm(freqHaploLi ~ Habitat*haplotype_Li, data = freqHaploLi)
car::Anova(LiLocusMod, type = 3)

# Model testing for variation in Ac haplotype relative frequency by haplotype and habitat type
LiLocusMod_Simp <- lm(simpson ~ Habitat, data = simpsDivLi)
summary(LiLocusMod_Simp, type = 3)

# Get Simpson's index in each habitat
simpsDivLi %>%
  group_by(Habitat) %>%
  summarize(meanSimp = mean(simpson))

#### TABLES ####

## TABLE 1 ##

# Linear cline models only. Part of table 1
linearClines <- linearClineModelOnly(city_dataframes)

# Extract number of populations and plants per city. Merge with linear clines model output above
table1 <- read.csv("data-clean/AllCities_AllPlants.csv") %>%
  group_by(City) %>%
  summarize(numPops = n_distinct(Population),
            numPlants = n()) %>%
  merge(., linearClines, by = "City")

# Write table 1 to disk
write_csv(table1, "analysis/tables/Table1_cityClineSummary.csv")

## TABLE 2 ##

table2 <- as.data.frame(rbind(
  c("Full"),
  summary(HCN_modAvg)$coefmat.full)) %>%
  rownames_to_column()

write_csv(table2, "analysis/tables/Table-2_freqHCN-FullModelAvg.csv")

# TABLE S2

tableS2 <- as.data.frame(rbind(
  c("Conditional"),
  summary(HCN_modAvg)$coefmat.subset)) %>%
  rownames_to_column()

write_csv(tableS2, "analysis/tables/Table-S2_freqHCN-CondModelAvg.csv")

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

## FIGURE 2 ##

colours = c("coral3", "cadetblue3", "burlywood3", "brown4",
            "blue1", "aquamarine2", "darkorchid3", "darkorange2",
            "khaki3", "hotpink2", "peru", "navy",
            "yellow3", "thistle3", "springgreen3", "cyan")
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

ggsave(filename = "analysis/figures/Figure2_HCN-by-distance.pdf", 
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

ggsave(filename = "analysis/figures/Figure3a_HCN-by-NumDaysNegNoSnow.pdf", 
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

ggsave(filename = "analysis/figures/Figure3b_HCN-by-PC1.pdf", 
       plot = HCN_by_PC1, device = 'pdf', units = 'in',
       width = 5, height = 5, dpi = 600)

# Figure 3b inset
envPCA_HCN_vars <- fviz_pca_var(envPCAHCN,
             labelsize = 6,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
) + ng1 + xlab("PC1 (90.2%)") + ylab("PC2 (7%)")

ggsave(filename = "analysis/figures/Figure3inset_envPCA_HCN_vars.pdf", 
       plot = envPCA_HCN_vars, device = 'pdf', units = 'in',
       width = 5, height = 5, dpi = 600)

## FIGURE 4 ##

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

ggsave(filename = "analysis/figures/Figure4_Slope-by-PC1.pdf", 
       plot = Slope_by_PC1, device = 'pdf', units = 'in',
       width = 5, height = 5, dpi = 600)

# Figure 3 inset
envPCA_Slope_vars <- fviz_pca_var(envPCAslope,
                                labelsize = 6,
                                col.var = "contrib", # Color by contributions to the PC
                                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                repel = TRUE     # Avoid text overlapping
) + ng1 + xlab("PC1 (92.8%)") + ylab("PC2 (3.7%)")

ggsave(filename = "analysis/figures/Figure4inset_envPCA_slope_vars.pdf", 
       plot = envPCA_Slope_vars, device = 'pdf', units = 'in',
       width = 5, height = 5, dpi = 600)

#### SUPPLEMENTARY FIGURES ####

## HAPLOTYPE FREQUENCIES

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

ggsave(filename = 'analysis/figures/supplemental/Li_haplotype_freqs.pdf', 
       plot = Li_Haplo_Freqs, device = "pdf", 
       width = 8, height = 5, dpi = 300)

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

ggsave(filename = 'analysis/figures/supplemental/Ac_haplotype_freqs.pdf', 
       plot = Ac_Haplo_Freqs, device = "pdf", 
       width = 8, height = 5, dpi = 300)

## INDIVIDUAL CLINE BIPLOTS

clineBiplot <- function(df, response_var, outpath, model_order_df){
  
  # Get city name
  city_name <- df$City[1]
  # print(city_name)
  if(model_order_df[model_order_df$City == city_name, response_var] == "linear"){
    
    plot <- df %>%
      ggplot(., aes_string(x = "std_distance", y = response_var)) +
      geom_point(colour = "black", size = 3.5) +
      geom_smooth(method = "lm", se = FALSE, colour = "black", size = 2) + 
      ylab("Frequency of HCN") + xlab("Standardized distance") +
      ng1
    
    # Full path to which data frame will be written
    path <- paste0(outpath, city_name, ".pdf")
    # print(path)
    
    # Write dataframe
    ggsave(filename = path, plot = plot, device = "pdf", 
           width = 5, height = 5, dpi = 300)
    
  }else if(model_order_df[model_order_df$City == city_name, response_var] == "quadratic"){
    
    plot <- df %>%
      ggplot(., aes_string(x = "std_distance", y = response_var)) +
      geom_point(colour = "black", size = 3.5) +
      geom_smooth(method = "lm", formula = y ~ x + I(x^2), 
                  se = FALSE, colour = "black", size = 2) + 
      ylab("Frequency of HCN") + xlab("Standardized distance") +
      ng1   
    
    path <- paste0(outpath, city_name, ".pdf")
    ggsave(filename = path, plot = plot, device = "pdf", 
           width = 5, height = 5, dpi = 300)
    
  }
}

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

## TEST FIGURES

haplotype_data %>%
  group_by(City, Habitat) %>%
  summarize(count = n())

#### CORRELATIONS AMONG WEATHER VARIABLES

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
write_csv(corr_mat, path = "analysis/weatherCorrMat.csv", col_names = TRUE)


  