
# Load required packages
library(dplyr)
library(broom)

# Load data with population-level data for all cities
datPops <- read.csv("data-clean/AllCities_AllPopulations.csv")

### FUNCTIONS

best_fit_model <- function(city, population_data, gene){
  
  # Define models based on whether frequency of HCN, Ac, or Li is being analyzed.
  if (gene == "HCN"){
    quadratic_model = lm(FreqC~Distance + Distance_squared, data = population_data) # Specify quadratic mode
    linear_model = update(quadratic_model, ~. -Distance_squared) # Specify linear model
  }else if (gene == "Ac"){
    quadratic_model = lm(FreqAc~Distance + Distance_squared, data = population_data) # Specify quadratic mode
    linear_model = update(quadratic_model, ~. -Distance_squared) # Specify linear model
  }else if (gene == "Li"){
    quadratic_model = lm(FreqLi~Distance + Distance_squared, data = population_data) # Specify quadratic mode
    linear_model = update(quadratic_model, ~. -Distance_squared) # Specify linear model
  }
  
  AIC_quad = AIC(quadratic_model) # Get AIC of quadratic model
  AIC_lin = AIC(linear_model) # Get AIC of linear model
  
  if (abs(AIC_quad) - abs(AIC_lin) > 2){ # If quadratic model AIC is > 2 from linear model AIC
    best_fit_model = quadratic_model # Then best fit model is quadratic
    best_fit = "quadratic"
  }else { # Otherwise (i.e. quadratic model is not better fit)
    best_fit_model = linear_model # Then best fit model is linear
    best_fit = "linear"
  }
  
  # Write best fit model output to disk
  output <- broom::tidy(best_fit_model)
  output_file <- sprintf("analysis/cline-model-output/%s-output/%s_cline-model_%s.csv", gene, city, gene)
  write.csv(output, output_file)
  
  # Return best fit model order
  return(best_fit)
}

### MODEL FIT AND CLINE MODEL OUTPUTS

# Here, I will assess whether the best fit model is linear or quadratic. Quadratic models are
# assumed to be a greater fit to the HCN, Ac or Li frequency data if they increase model
# AIC by more than 2 points. The best fit mode order (i.e. quadratic or linear) for each city and 
# gene (i.e. HCN, Ac, or Li) will be written to a dataframe. In addition, the best fit model 
# output for each city and gene will be written as a .csv file with all coefficients and P-values. 

# Get city names
Cities <- unique(datPops$City)

# Initialize dataset that will hold best fit model type
best_fit_model_dataset <- data.frame(City = character(),
                                    HCN = character(),
                                    Ac = character(),
                                    Li = character(),
                                    stringsAsFactors = FALSE)

# Iterate over cities
for (i in 1:length(Cities)){
  
  city = as.character(Cities[[i]]) # Extract city as character
  dat <- datPops %>% filter(City == city) # Subset dataset to to include only current city
  
  # These are the cities for which frequency data is available at the individual loci 
  # underlying HCN
  Gene_cities = c("Atlanta", "Baltimore", "Charlotte", 
                  "Cleveland", "Jacksonville", "NewYork", 
                  "Norfolk", "Toronto", "Washington")
  
  # Write best fit model output for HCN frequency to disk. Extract model order (quadratic or linear).
  best_fit_HCN <- best_fit_model(city, dat, "HCN")
  
  # If there is data for individual alleles
  if (city %in% Gene_cities){
    
    # Get model order of best fit model. Write output to disk
    best_fit_Ac <- best_fit_model(city, dat, "Ac") 
    best_fit_Li <- best_fit_model(city, dat, "Li")
    
    # Vector of model order for HCN, Ac, and Li to add to dataframe
    data_to_bind = c(city, best_fit_HCN, best_fit_Ac, best_fit_Li)
  }else{
    
    # "NA" for Ac and Li model order since underlying gene data unavailable
    data_to_bind = c(city, best_fit_HCN, "NA", "NA")
  }
  
  # Add model order vector to dataframe
  best_fit_model_dataset[i, ] = data_to_bind
}

write.csv(best_fit_model_dataset, "analysis/best_fit_model_order.csv")
