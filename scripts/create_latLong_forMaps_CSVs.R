# Load required packages
library(tidyverse)

datPops <- read_csv("data-clean/AllCities_AllPopulations.csv")
latLongs <- read_csv("data-raw/Lat-Longs_AllCities_AllPops.csv") %>%
  select(-Transect)
city_centres <- read_csv("data-raw/Lat-longs_City-centres.csv") %>%
  select(-Comments)

latLongs_forMap <- function(df, city_centres, latLongs, outpath){
  
  city_name <- df$City[1]

  centre <- city_centres %>% filter(City == city_name)
  df_out <- datPops %>%
    filter(City == city_name) %>%
    select(City, Population, freqHCN) %>%
    left_join(latLongs %>% filter(City == city_name)) %>%
    rbind(c(centre$City, "City Centre", "NA", centre$Latitude, centre$Longitude))
  
  path <- paste0(outpath, city_name, ".csv")
  write_csv(df_out, path)
}

# Create list with city dataframes as elements
city_df_list <- datPops %>% split(.$City)

outpath <- "GIS/city_maps/latLongs_CSVs/"
purrr::walk(city_df_list, latLongs_forMap,
            city_centres = city_centres, 
            latLongs = latLongs, outpath = outpath)


