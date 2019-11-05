# Script to perform power analyses to examine the minimum number of plants per population,
# and number of populations, to sample to identify urban-rural clines in HCN without loss of power. 


#####Load required packages#####
library(tidyverse)
library(grid)

set.seed(42)

# Load in data for Ken's cities
# Coltypes required to convert population to character
# due to presence of letters for some populations after first 1000 rows.
KT_plants <- read_csv("data-clean/AllCities_AllPlants.csv",
                      col_types = "ccciddiiicccd") %>% 
  filter(City %in% c("Toronto", "NewYork", "Boston", "Montreal"))

# Load in population level data
KT_pops <- read_csv("data-clean/AllCities_AllPopulations.csv") %>%
  filter(City %in% c("Toronto", "NewYork", "Boston", "Montreal"))

#Theme used for plots throughout script#
ng1=theme(aspect.ratio=0.7,panel.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border=element_blank(),
          axis.line.x = element_line(color="black",size=1), 
          axis.line.y = element_line(color="black",size=1),
          axis.ticks=element_line(color="black"), 
          axis.text=element_text(color="black",size=15), 
          axis.title=element_text(color="black",size=1), 
          axis.title.y=element_text(vjust=2,face="bold",size=17),
          axis.title.x=element_text(vjust=0.1,face="bold",size=17),
          axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          legend.position = "top", legend.direction="horizontal", 
          legend.text=element_text(size=17), legend.key = element_rect(fill = "white"), 
          legend.title = element_text(size=17),legend.key.size = unit(1.5, "cm"))

#### FUNCTIONS ####

# Sample plants and return slope
FunSlope <- function(dat, size){
  Plants_by_pop <- dat %>% 
    group_by(City, Population) %>% 
    # distinct(Plant) %>% 
    sample_n(size = size, replace = T) %>%
    select(City, Plant, HCN_Result, Distance, Population)
  
  datFreq <- Plants_by_pop %>%
    group_by(City, Population, Distance) %>%
    summarize(n = n(),
              n_hcn = sum(HCN_Result),
              freqHCN = n_hcn / n)
  
  model <- lm(freqHCN ~ Distance, data = datFreq)
  model_summary <- summary(model)
  slope <- model_summary$coefficients[2]
  return(slope)
}

# Sample populations, return slope
FunSlopePop<-function(dat, size){
  Pop_sample <- dat %>% group_by(City) %>% 
    sample_n(size = size, replace = T) %>%
    select(City, freqHCN, Distance, Transect, Population)
  
  
  model <- lm(freqHCN ~ Distance, data = Pop_sample)
  model_summary <- summary(model)
  slope <- model_summary$coefficients[2]
  return(slope)
}

calcProb <- function(num_vector, obs_val, nreps){
  
  prob <- (sum(abs(num_vector) >= (obs_val))) / nreps
  return(prob)
  
}

#### ANALYSIS OF BOSTON ####

# Model for Boston
KT_pops_Bos <- KT_pops %>% filter(City == "Boston")
ModBos <- lm(freqHCN ~ Distance, data = KT_pops_Bos)
SumBos <- summary(ModBos)
obs_BosSlope <- SumBos$coefficients[2]

# Get plant-level data for Boston, which shows the weakest cline overall
KT_plants_Bos <- KT_plants %>% filter(City == "Boston")

### SAMPLING PLANTS ###

# Slopes for Boston
S1B <- replicate(1000, FunSlope(KT_plants_Bos, size = 20))
S2B <- replicate(1000, FunSlope(KT_plants_Bos, size = 19))
S3B <- replicate(1000, FunSlope(KT_plants_Bos, size = 18))
S4B <- replicate(1000, FunSlope(KT_plants_Bos, size = 17))
S5B <- replicate(1000, FunSlope(KT_plants_Bos, size = 16))
S6B <- replicate(1000, FunSlope(KT_plants_Bos, size = 15))
S7B <- replicate(1000, FunSlope(KT_plants_Bos, size = 14))
S8B <- replicate(1000, FunSlope(KT_plants_Bos, size = 13))
S9B <- replicate(1000, FunSlope(KT_plants_Bos, size = 12))
S10B <- replicate(1000, FunSlope(KT_plants_Bos, size = 11))
S11B <- replicate(1000, FunSlope(KT_plants_Bos, size = 10))

#Bind results into table#
SlopeBosPlants <- cbind(S1B, S2B, S3B, S4B, S5B, S6B, S7B, S8B, S9B, S10B, S11B)
SlopeBosPlants <- as.data.frame(SlopeBosPlants)

plotBosSlopePlants <- SlopeBosPlants %>% 
  select_if(is.numeric) %>% 
  purrr::map_dfr(~data.frame(prob = calcProb(., obs_BosSlope, 1000))) %>% 
  mutate(num_plants = 20:10) %>% 
  ggplot(., aes(x = num_plants, y = prob)) +
  geom_point(size = 2, colour = "black") +
  coord_cartesian(ylim = c(0.9, 0.98)) + 
  scale_x_reverse(breaks = seq(from = 10, to = 20, by = 2)) +
  scale_y_continuous(breaks = seq(from = 0.9, to = 0.98, by = 0.02)) +
  geom_vline(xintercept = 15, linetype = "dashed") +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  xlab("Number of plants") + ylab("P(slope \u2265 observed Boston slope)") +
  ng1
plotBosSlopePlants

# Save plot to disk
ggsave("analysis/figures/supplemental/samplingPlants_slopeBoston.pdf",
       plot = plotBosSlopePlants, dpi = 300, width = 6, height = 6, units = "in",
       device = cairo_pdf)

