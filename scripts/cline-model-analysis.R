

### SETUP ###

# Clear environement, if necessary
rm(list=ls())
.rs.restartR()

# Load required packages
library(tidyverse)

# Load data with population-level data for all cities
datPops <- read.csv("data-clean/AllCities_AllPopulations.csv")

### FUNCTIONS ###

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

# Split cities into different dataframes stored as list
city_dataframes <- split(datPops, datPops$City)

# Write cline results to disk and global environment
writeClineResults(city_dataframes)

