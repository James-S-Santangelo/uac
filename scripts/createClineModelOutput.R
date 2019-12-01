# Predicting the strength of urban-rural clines in a 
# Mendelian polymorphism along a latitudinal gradient 
#
# Authors: James S. Santangelo, Ken A. Thompson, Beata Cohan
# Jibran Syed, Rob W. Ness, Marc T. J. Johnson
#
#
# Generate table with betas, p-values, t-stats, and std errors
# for best fit cline model for each city and for HCN, Ac, and Li.
# Also add the betas from a linear model only (i.e., regardless of fit)
# to each city. Ass logistic regression coefficients.

# Load in population data
datPops <- read_csv("data-clean/AllCities_AllPopulations.csv") %>% 
  mutate(Population = as.character(Population))

###################################################
#### STEP 1: SLOPES, ETC., FROM BEST FIT MODEL ####
###################################################

# For each city and response variable (i.e., HCN frequency,
# predicted Ac, or predicted Li frequency), we asses whether
# a first order or quadratic model fits the data best. We then 
# run this model and extract the beta coefficients for all terms,
# the standard error, the t-statistic, and the p-value. All of these
# values are renamed and concatenated into a single dataframe

# Possible values for the response variables
possibleResponses <- c("freqHCN", "AcHWE", "LiHWE")

mod_out <- c() # List to hold resulting dataframes

# Loop over response variables
for(res in possibleResponses){
  
  mod_out[[res]] <- datPops %>% 
    split(.$City) %>%  # Separate list for each city
    map_dfr(., ~tidy(runBestFitModel(res, .)), .id = "City") %>%  # Map function to get cline results for each city
    filter(term != "(Intercept)") %>%  
    gather(estimate:p.value, key = var, value = val) %>% 
    mutate(term = case_when(term == "std_distance" ~ "Lin",
                            term == "std_distance_squared" ~ "Quad"),
           var = case_when(var == "estimate" ~ "beta",
                           var == "std.error" ~ "stdErr",
                           var == "statistic" ~ "t",
                           var == "p.value" ~ "pval")) %>% 
    mutate(response = res) %>% 
    mutate(var = paste(paste0(var, term), response, sep = "_")) %>% 
    select(-term, -response) %>% 
    spread(key = var, value = val)
  
}

allModelOutputs <- reduce(mod_out, left_join, by = "City")

##################################################################
#### STEP 2: SLOPES FROM FIRST ORDER MODEL, REGARDLESS OF FIT ####
##################################################################

# Here we rund only a first order linear model for each city, 
# regardless of whether this provides the best fit to the data. 
# These are the slopes that should be compared between cities since
# they are measuring the same thing. 

# Run first order regression for each city
linSlopesOnly <- datPops %>% 
  group_by(City) %>% 
  do(modlm = lm(freqHCN ~ std_distance, data = .)) %>% 
  tidy(modlm, datPops) %>% 
  filter(term != "(Intercept)") %>% 
  select(City, estimate) %>% 
  rename("betaLinOnly" = estimate)


###################################################
#### STEP 3: LOGISTIC REGRESSION FOR EACH CITY ####
###################################################

# Here we perform a logistic regression for each city and
# extract the beta coefficients, standard error, z stat, 
# and P-value. 

# Load data for all plants 
datPlants <- read_csv("data-clean/AllCities_AllPlants.csv",
                      col_types = "ccncddnnnd") %>% 
  left_join(., datPops %>% select(City, Population, std_distance), by = c("City", "Population"))

# Logistic regression for each city
logRegs <- datPlants %>% 
  group_by(City) %>% 
  do(modlog = glm(HCN_Result ~ std_distance, family = "binomial", data = .)) %>% 
  tidy(modlog, datPlants) %>% 
  filter(term != "(Intercept)") %>% 
  select(-term) %>% 
  rename("betaLog" = estimate,
         "stdErrLog" = "std.error",
         "zLog" = "statistic",
         "pvalLog" = "p.value")

##########################################################
#### STEP 4: COMBINE DATAFRAMES WITH ALL COEFFICIENTS ####
##########################################################

allModelOutputs <- allModelOutputs %>% 
  left_join(., linSlopesOnly, by = "City") %>% 
  left_join(., logRegs, by = "City") %>% 
  mutate_if(is.numeric, round, 4)
write_csv(allModelOutputs, "analysis/clineModelOutput.csv")

