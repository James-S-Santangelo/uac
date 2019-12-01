# Predicting the strength of urban-rural clines in a 
# Mendelian polymorphism along a latitudinal gradient 
#
# Authors: James S. Santangelo, Ken A. Thompson, Beata Cohan
# Jibran Syed, Rob W. Ness, Marc T. J. Johnson
#
#
# Script cleans the raw Haplotype data for use in analyses.

# Load in haplotype data
haplo <- read_csv("data-raw/Haplotype-data.csv")

# Filter out columns and refactor levels
haplo_modified <- haplo %>%
  select(-aliquot_plate, -well_letter, -"well_#", -"sort column",
         -starts_with("Misprime"), -matches("comments|Comments")) %>%
  rename("HCN_Result" = "FA_HCN _identity") %>%
  separate(.,plant_sample, into=c("City","pop_plant"),sep = "(?<=[A-Z][A-Z])", remove = TRUE) %>%
  separate(., pop_plant, into = c("Population", "Plant"), sep = "[.]", remove = TRUE) %>%
  mutate(City = recode(City,
                       "AT" = "Atlanta",
                       "WA" = "Washington D.C.",
                       "BA" = "Baltimore",
                       "CL" = "Cleveland",
                       "NF" = "Norfolk",
                       "NY" = "NewYork",
                       "JA" = "Jacksonville")) %>%
  
  # Assign haplotypes for Ac locus
  mutate(haplotype_Ac = case_when(
    ac1.3 == 1 & ac2.2 == 1 & ac5.3 == 1 ~ "Ac", 
    ac1.3 == 0 & ac2.2 == 1 & ac5.3 == 1 ~ "ac2.2",
    ac1.3 == 0 & ac2.2 == 0 & ac5.3 == 1 ~ "ac5.3",
    # ac1.3 == 0 & ac2.2 == 1 & ac5.3 == 0 ~ "dc",
    TRUE ~ "unk"
  )) %>%
  
  # Assign haplotypes for Li locus
  mutate(haplotype_Li = case_when(
    li4.0 == 1 & li7.0 == 1 & li9.3 == 1 ~ "Li",
    li4.0 == 0 & li7.0 == 1 & li9.3 == 1 ~ "li6.6", 
    li4.0 == 0 & li7.0 == 1 & li9.3 == 0 ~ "li6.6DD", 
    li4.0 == 0 & li7.0 == 0 & li9.3 == 1 ~ "li9.1", 
    li4.0 == 0 & li7.0 == 0 & li9.3 == 0 ~ "li>11", 
    TRUE ~ "unk"
  )) %>%
  
  # Create column for overall plant haplotype
  mutate(haplotype_all = paste(haplotype_Ac, haplotype_Li, sep = "-"))

# Load dataset to retrieve distances
population_distances <- read_csv("data-clean/AllCities_AllPopulations.csv") %>%
  select(City, Population, Distance)

# Add distances and habitat type to haplotype data
haplo_modified <- haplo_modified %>%
  left_join(., population_distances, by = c("City", "Population")) %>%
  group_by(City) %>%
  mutate(Habitat = ifelse(Distance < 20, "Urban", "Rural"))

# Write clean haplotype data to disk
write_csv(haplo_modified, "data-clean/haplotypeData.csv")