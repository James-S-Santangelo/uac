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
#' @param df Dataset to use for running cline model.
#' @return Model order ("linear" or "quadratic") for best fit cline
getBestFitClineModelOrder <- function(response, df) {

  # Pull variables for linear model from dataset
  response_var <- df %>% pull(response) #[, response]
  std_distance <- df %>% pull("std_distance")
  std_distance_squared <- df %>% pull("std_distance_squared")

  if(all(is.na(response_var)) == TRUE){
    order <- "NA"
  }else{

    # Define linear and quadratic models for analysis of HCN with distance
    quadratic_model <-  lm(response_var ~ std_distance + std_distance_squared) # Specify quadratic model
    linear_model <-  update(quadratic_model, ~ . - std_distance_squared) # Specify linear model

    AIC_quad <- AIC(quadratic_model) # Get AIC of quadratic model
    AIC_lin <- AIC(linear_model) # Get AIC of linear model

    if (abs(AIC_quad) - abs(AIC_lin) > 2) {
      # If quadratic model AIC is > 2 from linear model AIC
      # bestFitClineModel <- quadratic_model # Then best fit model is quadratic
      order <- "quadratic"
    } else {
      # Otherwise (i.e. quadratic model is not better fit)
      # bestFitClineModel <- linear_model # Then best fit model is linear
      order <- "linear"
    }
  }


  # Return list with best fit model and model order as string
  return(order)
}


#' Runs model based on order (linear or quadratic) that fits best
#'
#' Quadratic model is assumed to be better fit if it improves model AIC
#' by more than 2 points.
#'
#' @param response Response variable to use in cline models. One of "freqHCN",
#' "AcHWE", or "LiHWE".
#' @param df Dataset to use for running cline model.
#' @return 'lm' model object for best fit model
runBestFitModel <- function(response, df){

  city <- as.character(unique(df$City))

  # Pull variables for linear model from dataset
  response_var <- df %>% pull(response)
  std_distance <- df %>% pull("std_distance")
  std_distance_squared <- df %>% pull("std_distance_squared")

  # order <- getBestFitClineModelOrder(response, df) # Get model order

  if(all(is.na(response_var)) == TRUE){
    model = "NA"
  }else{
    order <- getBestFitClineModelOrder(response, df)
    if(order == "linear"){
      model <- lm(response_var ~ std_distance) # Only first order term
      write_csv(tidy(model),
                path = sprintf("analysis/inividual-cline-models/%s/%s-%s-cline-model.csv",
                               response, city, response),
                append = FALSE)
    }else if(order == "quadratic"){
      model <- lm(response_var ~ std_distance + std_distance_squared) # First and second order terms
      write_csv(tidy(model),
                path = sprintf("analysis/inividual-cline-models/%s/%s-%s-cline-model.csv",
                               response, city, response),
                append = FALSE)
    }
    return(model)
  }
}

#' Generates biplot of frequency of HCN, Ac, or Li, against distance
#'
#' @param df Dataframe with population mean allele or phenotype frequencies and distance
#' @param response_var One of 'freqHCN', 'AcHWE', or 'LiHWE'
#' @param outpath Path to which figures should be written
#' @param model_order_df Dataframe specifying order of best fit linear model as
#'    'linear' or 'quadratic'. Gerenated with `getBestFitClineModel`
#'
#' @return None. Saves biplot to disk
clineBiplot <- function(df, response_var, outpath, model_order_df){

  # Get city name
  city_name <- as.character(unique(df$City))

  # Get model order
  order <- getBestFitClineModelOrder(response_var, df)

  if(order == "linear"){

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

  }else if(order == "quadratic"){

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


#' Plot logistic regression plot for each city
#' 
#' @param df_allPlants Dataframe with individual plant phenotype data
#' @param city City for which logistic regression plot should be created
#' @param tag Tag specifying position in multipannelled figure (e.g., A, B, D, or D)
#' 
#' @return None. Saves logistic regression plot to disk.
plotLogReg <- function(df_allPlants, city, tag){

  df <- df_allPlants %>% filter(City == city)

  mod <- glm(HCN_Result ~ std_distance, family = "binomial", data = df)
  beta_val <- round(summary(mod)$coefficients[2], 3)
  pval <- summary(mod)$coefficients[8]
  pval <- ifelse(pval < 0.001, "< 0.001", round(pval, 3))
  # print(c(coef, pval))

  plot <- ggplot(df, aes(x=std_distance, y=HCN_Result)) +
    geom_point(alpha=.5) +
    stat_smooth(method="glm", se=TRUE, fullrange=TRUE,
                method.args = list(family=binomial),
                color = "black") +
    ylab("Presence/absence of HCN") +
    xlab("Standardized distance 
to urban center") +
    labs(tag = tag) +
    coord_cartesian(xlim = c(0, 1.05)) +
    scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.25)) +
    annotate(geom = "text", label = city,
             x = (max(df$std_distance) - min(df$std_distance)) / 2,
             y = ifelse(city == "Jacksonville", 0.75, 0.85),
             size = 5) +
    annotate(geom = "text", parse = TRUE,
             label = paste("italic(beta)==", beta_val),
             x = ifelse(city == "Jacksonville", 1, 0),
             y = ifelse(city == "Jacksonville", 0.5, 0.85),
             size = 4,
             hjust = ifelse(city == "Jacksonville", "right", "left")) +
    annotate(geom = "text", parse = TRUE,
             label = ifelse(is.character(pval),
                            paste("italic(P)", pval),
                            paste("italic(P)==", pval)),
             x = ifelse(city == "Jacksonville", 1, 0),
             y = ifelse(city == "Jacksonville", 0.4, 0.75),
             size = 4,
             hjust = ifelse(city == "Jacksonville", "right", "left")) +
    ng1 + theme(plot.tag = element_text(size = 15),
                plot.tag.position = c(0.15, 0.83),
                plot.margin = margin(-2, 0.5, -2, 0, "cm"),
                axis.title.y = element_text(vjust = 2, size = 15),
                axis.title.x = element_text(vjust = 0.1, size = 15))

  return(plot)
}

#' Convert degrees to radians
#' 
#' @param deg Degrees
#' 
#' @return Radians
deg2rad <- function(deg) return(deg*pi/180)

#' Calculates the geodesic distance between two points 
#' 
#' @description Geodesic distance between two points is specified by radian
#'   latitude/longitude using the Haversine formula (hf)
#'   
#' @param long1 Longitude of first coordinate in radians
#' @param lat1 Latitude of first coordinate in radians
#' @param long2 Longitude of second coordinate in radians
#' @param lat2 Latitude of second coordinate in radians
#' 
#' @return Distance between both coordinates in kilometers
haversine <- function(long1, lat1, long2, lat2) {

  # Ensure Lats and Longs are in radians
  long1 <- deg2rad(long1)
  lat1 <- deg2rad(lat1)
  long2 <- deg2rad(long2)
  lat2 <- deg2rad(lat2)

  # Calculate geodesic distance based on havesine formala
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  d = R * c
  return(d) # Distance in km
}

#' Creates simplified lat/long dataframe for creating maps in QGIS
#' 
#' @param df Dataframe with population means
#' @param city_centres Dataframe with the coordinates of city centres for each city
#' @param latLongs Dataframe with lat longs for each populations
#' 
#' @return Dataframe with lat longs of populations and city centres for each city.
latLongs_forMap <- function(df, city_centres, latLongs){

  city_name <- df$City[1]

  centre <- city_centres %>% filter(City == city_name)
  df_out <- datPops %>%
    filter(City == city_name) %>%
    select(City, Transect, Population, freqHCN) %>%
    left_join(., latLongs %>% filter(City == city_name), by = c("City", "Population")) %>%
    rbind(c(centre$City, "NA","City Centre", "NA", centre$Latitude, centre$Longitude))

  return(df_out)
}


#' Randomly sample plants within populations
#' 
#' @description Randomly sample plants within populations, run linear moden of
#'   mean HCN frequency against distance, and return the slope of this model
#'   
#' @param dat Population mean dataframe
#' @param size Number of plants to sample
#' 
#' @return Slope of linear model
FunSlope <- function(dat, size){
  Plants_by_pop <- dat %>% 
    group_by(City, Population) %>% 
    # distinct(Plant) %>% 
    sample_n(size = size, replace = T) %>%
    select(City, Plant, HCN_Result, Distance, Population)
  
  datFreq <- Plants_by_pop %>%
    group_by(City, Population, Distance) %>%
    summarise(n = n(),
              n_hcn = sum(HCN_Result),
              freqHCN = n_hcn / n)
  
  model <- lm(freqHCN ~ Distance, data = datFreq)
  model_summary <- summary(model)
  slope <- model_summary$coefficients[2]
  return(slope)
}

#' Calculate P-value from null distribution of slopes
#' 
#' @description Calculates probability that observed value is greater or equal
#'   to expected value based on randomly generated null distribution
#'   
#' @param num_vector Vector with test statistics generated by randomly
#'   reshuffling data
#' @param obs_val Observed value of test statistic
#' @param nreps Number of resamplings of dataframe
#' 
#' @return P-value 
calcProb <- function(num_vector, obs_val, nreps){
  
  prob <- (sum(abs(num_vector) >= (obs_val))) / nreps
  return(prob)
  
}
