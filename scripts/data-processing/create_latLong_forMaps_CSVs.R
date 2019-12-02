# Predicting the strength of urban-rural clines in a
# Mendelian polymorphism along a latitudinal gradient
#
# Authors: James S. Santangelo, Ken A. Thompson, Beata Cohan
# Jibran Syed, Rob W. Ness, Marc T. J. Johnson
#
#
# Script to generate dataset with Latitudes and Longitudes for each city and population
# with city centres. Used for creating maps.

# Load datasets
datPops <- read_csv("data-clean/AllCities_AllPopulations.csv") %>%
  mutate(Population = as.character(Population))
latLongs <- read_csv("data-raw/Lat-Longs_AllCities_AllPops.csv") %>%
  select(-Transect) %>%
  mutate(Population = as.character(case_when(City =="Boston" & Population == "4A" ~ "44",
                                             City =="Boston" & Population == "4B" ~ "45",
                                             City =="NewYork" & Population == "11-I" ~ "11",
                                             City =="NewYork" & Population == "11-II" ~ "12",
                                             TRUE ~ Population)))
city_centres <- read_csv("data-raw/Lat-longs_City-centres.csv") %>%
  select(-Comments)

# Create list with city dataframes as elements
city_df_list <- datPops %>% split(.$City)

# Generate single dataframe with populations and City centres
allPops_withCentres <- purrr::map_dfr(city_df_list, latLongs_forMap,
               city_centres = city_centres,
               latLongs = latLongs)
write_csv(allPops_withCentres, "data-clean/allPops_withCentres.csv")

