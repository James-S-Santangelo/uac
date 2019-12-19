# Predicting the strength of urban-rural clines in a 
# Mendelian polymorphism along a latitudinal gradient 
#
# Authors: James S. Santangelo, Ken A. Thompson, Beata Cohan
# Jibran Syed, Rob W. Ness, Marc T. J. Johnson
#
#
# Script to carry out all analyses reported in manuscript

###############
#### SETUP ####
###############

# Change default contrasts to enable type III SS
options(contrasts = c("contr.sum", "contr.sum"))
# options(contrasts = c("contr.treatment", "contr.poly")) # Default

# Clear environement, if necessary
# rm(list=ls())
# .rs.restartR()

## LOAD DATASETS FOR ANALYSES ##

# Population-mean HCN frequencies
datPops <- read_csv("data-clean/AllCities_AllPopulations.csv") %>% 
  mutate(Population = as.character(Population))

# Individual plant-level data. Append std_distance
datPlants <- read_csv("data-clean/AllCities_AllPlants.csv",
                      col_types = "ccncddnnnnd") %>% 
  left_join(., datPops %>%select(City, Population, std_distance), by = c("City", "Population"))

# Data with environmental and model summaries by city
citySummaryData <- read_csv("data-clean/citySummaryData.csv") 

# Load in data
haplotype_data <- read_csv("data-clean/haplotypeData.csv")


##############################################
#### ANALYSIS OF CLINES ACROSS ALL CITIES ####
##############################################

## POPULATION-MEAN HCN FREQUENCIES ###

# Anova for overall effect of distance, city, and whether the strenth of clines
# varies across cities (i.e., distance x city interaction)
clinesAllCities <- lm(freqHCN ~ std_distance*City, data = datPops)
summary(clinesAllCities)
Anova(clinesAllCities, type = 3)

## INDIVIDUAL PLANT PHENOTYPE DATA ##

# Perform logistic regression. P-values from type III SS. 
logClinesAllCitiesFreqs <- glm(freqHCN ~ std_distance*City, family = "binomial", data = datPops, weights = n_HCN)
summary(logClinesAllCitiesFreqs)
Anova(logClinesAllCitiesFreqs, type = 3)

# To get predicted probabilities, averaged across all cities
logClinesAllCitiesFreqs_distOnly <- glm(freqHCN ~ std_distance, family = "binomial", data = datPops, weights = n_HCN)
intercept <- coef(logClinesAllCitiesFreqs_distOnly)["(Intercept)"]
distance_effect <- coef(logClinesAllCitiesFreqs_distOnly)["std_distance"]

# Predicted probabilities.
probCyan_urban <- exp(intercept + distance_effect * 0) / (1 + (exp(intercept + distance_effect * 0)))
probCyan_rural <- exp(intercept + distance_effect * 1) / (1 + (exp(intercept + distance_effect * 1)))
(probCyan_rural - probCyan_urban) / probCyan_urban # Relative increase in probability of being cyanogenic from rural --> urban


###############################################
#### CORRELATIONS AMONG CLIMATIC VARIABLES ####
###############################################

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
  rename(PC1_SlopeLin = Dim.1) %>%
  rownames_to_column("City") %>%
  merge(., citySummaryDataForAnalysis, by = "City")

SlopeModLin <- lm(betaLinOnly ~ PC1_SlopeLin, data = citySummaryDataForAnalysis)
summary(SlopeModLin)

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
  rename(PC1_SlopeLog = Dim.1) %>%
  rownames_to_column("City") %>%
  merge(., citySummaryDataForAnalysis, by = "City")

SlopeModLog <- lm(betaLog ~ PC1_SlopeLog, data = citySummaryDataForAnalysis)
summary(SlopeModLog)

################################################
#### ANALYSIS OF HAPLOTYPES BY HABITAT TYPE ####
################################################

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
AcLocusMod_Rich <- lm(richness ~ Habitat, data = richnessHaplotAc)
summary(AcLocusMod_Rich)

# Mean
richnessHaplotAc %>% group_by(Habitat) %>% summarise(mean = mean(richness),
                                                     se = sd(richness) / sqrt(n()))

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

# Model testing for variation in "li" haplotype richness by habitat
LiLocusMod_Rich <- lm(richness ~ Habitat, data = richnessHaplotLi)
summary(LiLocusMod_Rich)

# Mean
richnessHaplotLi %>% group_by(Habitat) %>% summarise(mean = mean(richness),
                                                     se = sd(richness) / sqrt(n()))

###########################
#### HERBIVORY ANLYSES ####
###########################

# Analysis of changes in population-mean herbivory with distance for 4 cities
herbMods <- datPlants %>% 
  filter(!(is.na(Dmg.Avg))) %>% 
  group_by(City, Population, Distance, std_distance) %>% 
  summarise(meanHerb = mean(Dmg.Avg)) %>% 
  ungroup() %>% 
  group_by(City) %>% 
  do(mod = lm(meanHerb ~ Distance, data = .)) %>% 
  tidy(., mod) %>% 
  filter(term != "(Intercept)") %>% 
  mutate_if(is.numeric, round, 3)

# Analysis of changes in Average herbivory among only acyanogenic plants
herbMods2 <- datPlants %>% 
  filter(!(is.na(Dmg.Avg)) &
           HCN_Result == 0) %>% 
  group_by(City) %>% 
  do(mod = lm(Dmg.Avg ~ Distance, data = .)) %>% 
  tidy(., mod) %>% 
  filter(term != "(Intercept)") %>% 
  mutate_if(is.numeric, round, 3)
