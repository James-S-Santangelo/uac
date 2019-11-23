# Predicting the strength of urban-rural clines in a 
# Mendelian polymorphism along a latitudinal gradient 
#
# Authors: James S. Santangelo, Ken A. Thompson, Beata Cohan
# Jibran Syed, Rob W. Ness, Marc T. J. Johnson
#
#
# Functions used during analyses and figure/table generation


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
getBestFitClineModel <- function(response, data) {
  
  # Pull variables for linear model from dataset
  response <- data %>% pull(response) #[, response]
  std_distance <- data %>% pull("std_distance")
  std_distance_squared <- data %>% pull("std_distance_squared")
  
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

#' Retrieves the slopes and P-values for linear and quadratic (if present) terms from lm object
#'
#' @param best_model_output lm object representing best fit cline model
#' @param best_model_order String representing order of best fit model
#' @param clines_results Vector storing coefficients from model output
#' @return clines_results: Vector storing coefficients from model output
getClineModelOutputFromOrder <-
  function(best_model_output,
           best_model_order,
           clines_results) {
    
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

#' Writes cline model output for HCN, Ac, and Li to table
#'
#' @param dataframe_list A list containing individual dataframes with
#' the frequency of HCN, Ac (if present), and Li (if present), and
#' standardized distance to the urban centre and distance squared
#' (for quadratic cline model)
#' @return None: Writes dataframe to disk and assigns to global environment
writeClineResults <- function(dataframe_list) {
  
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
            "analysis/tables/supplemental/TableS3_clineOrderData.csv",
            na = "NA",
            append = FALSE)
}

#' Writes cline model output of linear model (first-order only)for HCN, Ac, and Li to table
#'
#' @param dataframe_list A list containing individual dataframes with
#' the frequency of HCN, Ac (if present), and Li (if present), and
#' standardized distance to the urban centre and distance squared
#' (for quadratic cline model)
#' @return modelOutputData A dataframe with slopes and P-values of first-order linear models
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

#' Generates biplot of frequency of HCN, Ac, or Li, against distance
#'
#' @param df Dataframe with population mean allele or phenotype frequencies and distance
#' @param response_var One of 'freqHCN', 'AcHWE', or 'LiHWE'
#' @param outpath Path to which figures should be written
#' @param model_order_df Dataframe specifying order of best fit linear model as
#'    'linear' or 'quadratic'. Gerenated with `getBestFitClineModel`
#' 
#' @return modelOutputData A dataframe with slopes and P-values of first-order linear models
clineBiplot <- function(df, response_var, outpath, model_order_df){
  
  # Get city name
  city_name <- df$City[1]
  # print(city_name)
  if(model_order_df[model_order_df$City == city_name, response_var] == "linear"){
    
    plot <- df %>%
      ggplot(., aes_string(x = "std_distance", y = response_var)) +
      geom_point(colour = "black", size = 3.5) +
      geom_smooth(method = "lm", se = FALSE, colour = "black", size = 2) + 
      ylab(sprintf("Frequency of %s", response_var)) + xlab("Standardized distance") +
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
      ylab(sprintf("Frequency of %s", response_var)) + xlab("Standardized distance") +
      ng1   
    
    path <- paste0(outpath, city_name, ".pdf")
    ggsave(filename = path, plot = plot, device = "pdf", 
           width = 5, height = 5, dpi = 300)
    
  }
}


plotLogReg <- function(df_allPlants, city, tag){
  
  df <- df_allPlants %>% filter(City == city)
  
  mod <- glm(HCN_Result ~ std_distance, family = "binomial", data = df)
  beta_val <- round(summary(mod)$coefficients[2], 3)
  pval <- summary(mod)$coefficients[8] 
  pval <- ifelse(pval < 0.001, "< 0.001", round(pval, 3))
  # print(c(coef, pval))
  
  plot <- ggplot(df, aes(x=Distance, y=HCN_Result)) + 
    geom_point(alpha=.5) +
    stat_smooth(method="glm", se=TRUE, fullrange=TRUE, 
                method.args = list(family=binomial),
                color = "black") + 
    ylab("") +
    xlab("") + 
    labs(tag = tag) + 
    annotate(geom = "text", label = city, 
             x = (max(df$Distance) - min(df$Distance)) / 2, 
             y = ifelse(city == "Jacksonville", 0.75, 0.85),
             size = 5) + 
    annotate(geom = "text", parse = TRUE, 
             label = paste("italic(beta)==", beta_val),
             x = ifelse(city == "Jacksonville", 40, 0),
             y = ifelse(city == "Jacksonville", 0.5, 0.85),
             size = 4,
             hjust = ifelse(city == "Jacksonville", "right", "left")) +
    annotate(geom = "text", parse = TRUE, 
             label = ifelse(is.character(pval), 
                            paste("italic(P)", pval),
                            paste("italic(P)==", pval)),
             x = ifelse(city == "Jacksonville", 40, 0), 
             y = ifelse(city == "Jacksonville", 0.4, 0.75),
             size = 4,
             hjust = ifelse(city == "Jacksonville", "right", "left")) +
    ng1 + theme(plot.tag = element_text(size = 15),
                plot.tag.position = c(0.15, 0.83),
                plot.margin = margin(-2, 0, -2, 0, "cm"))
  
  return(plot)
}
