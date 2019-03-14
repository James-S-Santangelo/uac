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

# Sample plants and return P-value
FunPval <- function(dat, size){
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
  p_val <- model_summary$coefficients[8]
  return(p_val)
}

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

# Sample populations, return P-value
FunPvalPop <- function(dat, size){
  Pop_sample <- dat %>% group_by(City) %>% 
    sample_n(size = size, replace = T) %>%
    select(City, freqHCN, Distance, Transect, Population)
  
  
  model <- lm(freqHCN ~ Distance, data = Pop_sample)
  model_summary <- summary(model)
  p_val <- model_summary$coefficients[8]
  return(p_val)
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

#### ANALYSIS OF BOSTON ####

# Model for Boston
KT_pops_Bos <- KT_pops %>% filter(City == "Boston")
ModBos <- lm(freqHCN ~ Distance, data = KT_pops_Bos)
SumBos <- summary(ModBos)

# Get plant-level data for Boston, which shows the weakest cline overall
KT_plants_Bos <- KT_plants %>% filter(City == "Boston")

### SAMPLING PLANTS ###

# Pvals
P1B <- replicate(1000, FunPval(KT_plants_Bos, size = 20))    
P2B <- replicate(1000, FunPval(KT_plants_Bos, size = 19))     
P3B <- replicate(1000, FunPval(KT_plants_Bos, size = 18))    
P4B <- replicate(1000, FunPval(KT_plants_Bos, size = 17))
P5B <- replicate(1000, FunPval(KT_plants_Bos, size = 16))
P6B <- replicate(1000, FunPval(KT_plants_Bos, size = 15))
P7B <- replicate(1000, FunPval(KT_plants_Bos, size = 14))
P8B <- replicate(1000, FunPval(KT_plants_Bos, size = 13))
P9B <- replicate(1000, FunPval(KT_plants_Bos, size = 12))
P10B <- replicate(1000, FunPval(KT_plants_Bos, size = 11))
P11B <- replicate(1000, FunPval(KT_plants_Bos, size = 10))

#Bind results into table#
PvalBosPlants <- cbind(P1B, P2B, P3B, P4B, P5B, P6B, P7B, P8B, P9B, P10B, P11B)
PvalBosPlants <- as.data.frame(PvalBosPlants)
Plants <- c(20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10)
MeansPvalBosPlants <- colMeans(PvalBosPlants)
BosPlantsPlotPval <- as.data.frame(cbind(Plants, MeansPvalBosPlants))

#Plot number of plants vs. mean P-value for Boston#
BosPlantsPlotPval <- ggplot(BosPlantsPlotPval, aes(x = Plants, y = MeansPvalBosPlants)) +
  geom_line(colour = "black",size = 1) + ylab("Mean P-value") + 
  xlab("Number of plants")  +geom_point(size = 2)+
  coord_cartesian(ylim = c(0.155, 0.24)) +
  scale_x_continuous(breaks = seq(from = 10, to = 20, by = 2)) +
  scale_y_continuous(breaks = seq(from = 0.16, to = 0.24, by = 0.02)) +
  ng1
BosPlantsPlotPval

# Save plot to disk
ggsave("analysis/figures/supplemental/samplingPlants_pvalBoston.pdf",
       plot = BosPlantsPlotPval, dpi = 300, width = 5, height = 5, units = "in")

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
Plants <- c(20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10)
MeansBosSlopePlants <- colMeans(SlopeBosPlants)
SdBosSlopePlants <- sapply(SlopeBosPlants, function(cl) sds = sd(cl))
BosPlantsPlotSlope <- as.data.frame(cbind(Plants, SdBosSlopePlants, MeansBosSlopePlants))

#Plot number of plants vs. mean slope for Boston#
plotBosSlopePlants <- ggplot(BosPlantsPlotSlope, aes(x = Plants, y = SdBosSlopePlants)) +
  geom_line(colour = "black",size = 1) + ylab("Mean sd of slope") + 
  xlab("Number of plants") + geom_point(size = 2) + 
  coord_cartesian(ylim = c(0.0009, 0.0014))+
  scale_y_continuous(breaks = seq(from = 0.0009, to = 0.0014, by = 0.0001))+
  scale_x_continuous(breaks = seq(from = 10, to = 20, by = 2)) +
  ng1
plotBosSlopePlants 

# Save plot to disk
ggsave("analysis/figures/supplemental/samplingPlants_slopeBoston.pdf",
       plot = plotBosSlopePlants, dpi = 300, width = 5, height = 5, units = "in")

### SAMPLING POPULATIONS ###

datBosFreq <- KT_pops %>% filter(City == "Boston")

#PVALS
P1B <- replicate(1000, FunPvalPop(datBosFreq, size = 50))    
P2B <- replicate(1000, FunPvalPop(datBosFreq, size = 45))     
P3B <- replicate(1000, FunPvalPop(datBosFreq, size = 40))    
P4B <- replicate(1000, FunPvalPop(datBosFreq, size = 35))
P5B <- replicate(1000, FunPvalPop(datBosFreq, size = 30))
P6B <- replicate(1000, FunPvalPop(datBosFreq, size = 25))
P7B <- replicate(1000, FunPvalPop(datBosFreq, size = 20))
P8B <- replicate(1000, FunPvalPop(datBosFreq, size = 15))
P9B <- replicate(1000, FunPvalPop(datBosFreq, size = 10))

#Bind results into table#
PvalBosPop <- cbind(P1B, P2B, P3B, P4B, P5B, P6B, P7B, P8B, P9B)
PvalBosPop <- as.data.frame(PvalBosPop)
Pops <- c(50, 45, 40, 35, 30, 25, 20, 15, 10)
MeansPvalBosPop <- colMeans(PvalBosPop)
SdPvalBosPop <- sapply(PvalBosPop, function(cl) sds = sd(cl))
df_BosPlotPvalPop <- as.data.frame(cbind(Pops, MeansPvalBosPop, SdPvalBosPop))

#Plot number of plants vs. mean P-value for Boston#
BosPlotPvalPop <- ggplot(df_BosPlotPvalPop, aes(x = Pops, y = MeansPvalBosPop)) +
  geom_line(colour = "black", size = 1) + ylab("Mean P-value") + 
  xlab("Number of populations") + geom_point(size = 2) + 
  coord_cartesian(ylim = c(0.1, 0.4)) +
  scale_y_continuous(breaks = seq(from = 0.1, to = 0.4, by = 0.05))+
  scale_x_continuous(breaks = seq(from = 10, to = 50, by = 5)) +
  ng1
BosPlotPvalPop 

# Save plot to disk
ggsave("analysis/figures/supplemental/samplingPops_pvalBoston.pdf",
       plot = BosPlotPvalPop, dpi = 300, width = 5, height = 5, units = "in")

# Slope
S1B <- replicate(1000, FunSlopePop(datBosFreq, size = 50))
S2B <- replicate(1000, FunSlopePop(datBosFreq, size = 45))
S3B <- replicate(1000, FunSlopePop(datBosFreq, size = 40))
S4B <- replicate(1000, FunSlopePop(datBosFreq, size = 35))
S5B <- replicate(1000, FunSlopePop(datBosFreq, size = 30))
S6B <- replicate(1000, FunSlopePop(datBosFreq, size = 25))
S7B <- replicate(1000, FunSlopePop(datBosFreq, size = 20))
S8B <- replicate(1000, FunSlopePop(datBosFreq, size = 15))
S9B <- replicate(1000, FunSlopePop(datBosFreq, size = 10))

#Bind results into table#
SlopeBosPop <- cbind(S1B, S2B, S3B, S4B, S5B, S6B, S7B, S8B, S9B)
SlopeBosPop <- as.data.frame(SlopeBosPop)
Pops <- c(50, 45, 40, 35, 30, 25, 20, 15, 10)
MeansBosSlopePop <- colMeans(SlopeBosPop)
SdBosSlopePop <- sapply(SlopeBosPop, function(cl) sds = sd(cl))
df_BosPlotSlopePop <- as.data.frame(cbind(SdBosSlopePop, Pops, MeansBosSlopePop))

#Plot Sd slope vs. number of populations#
plotBosSdSlopePop <- ggplot(df_BosPlotSlopePop, aes(x = Pops, y = SdBosSlopePop)) +
  geom_line(colour = "black",size = 1) + ylab("Mean sd of slope") + 
  xlab("Number of populations") + geom_point(size = 2) + 
  coord_cartesian(ylim = c(0.0010,0.0035)) +
  scale_y_continuous(breaks = seq(from = 0.001, to = 0.0035, by = 0.00035))+
  scale_x_continuous(breaks = seq(from = 10, to = 50, by = 5)) + 
  ng1
plotBosSdSlopePop 

# Save plot to disk
ggsave("analysis/figures/supplemental/samplingPops_slopeBoston.pdf",
       plot = plotBosSdSlopePop, dpi = 300, width = 5, height = 5, units = "in")

#### ANALYSIS OF TORONTO ####

# Use population data to figure out which Toronto transect showed strongest cline
KT_pops_TorA <- KT_pops %>% filter(City == "Toronto" & Transect == "A")
ModTorA <- lm(freqHCN ~ Distance, data = KT_pops_TorA)
SumTorA <- summary(ModTorA)
SumTorA$coefficients[8]

KT_pops_TorB <- KT_pops %>% filter(City == "Toronto" & Transect == "B")
ModTorB <- lm(freqHCN ~ Distance, data = KT_pops_TorB)
SumTorB <- summary(ModTorB)
SumTorB$coefficients[8]

KT_pops_TorC <- KT_pops %>% filter(City == "Toronto" & Transect == "C")
ModTorC <- lm(freqHCN ~ Distance, data = KT_pops_TorC)
SumTorC <- summary(ModTorC)
SumTorC$coefficients[8]


### SAMPLING PLANTS ###

#The strongest pattern is found in transect B of Toronto. Will use this for analysis#
# Get plant-level data for Toronto plants along transect B
KT_plants_TorB <- KT_plants %>% filter(City == "Toronto" & Transect == "B")

# Sample different numbers of plants with replacement from within populations
P1T <- replicate(1000, FunPval(KT_plants_TorB, size = 20))    #20 plants per population, with replacement
P2T <- replicate(1000, FunPval(KT_plants_TorB, size = 19))    #19 of plants per population, with replacement
P3T <- replicate(1000, FunPval(KT_plants_TorB, size = 18))    #18 of plants per population, with replacement
P4T <- replicate(1000, FunPval(KT_plants_TorB, size = 17))    #Etc...
P5T <- replicate(1000, FunPval(KT_plants_TorB, size = 16))
P6T <- replicate(1000, FunPval(KT_plants_TorB, size = 15))
P7T <- replicate(1000, FunPval(KT_plants_TorB, size = 14))
P8T <- replicate(1000, FunPval(KT_plants_TorB, size = 13))
P9T <- replicate(1000, FunPval(KT_plants_TorB, size = 12))
P10T <- replicate(1000, FunPval(KT_plants_TorB, size = 11))
P11T <- replicate(1000, FunPval(KT_plants_TorB, size = 10))

#Bind results into table#
PvalTor <- cbind(P1T, P2T, P3T, P4T, P5T, P6T, P7T, P8T, P9T, P10T, P11T)
PvalTor <- as.data.frame(PvalTor)
Plants <- c(20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10)
MeansTorPval <- colMeans(PvalTor)
TorPlotPval <- as.data.frame(cbind(Plants, MeansTorPval))

#Plot Means P-value vs. nnumber of plants for Toronto#
plotTorPVal <- ggplot(TorPlotPval, aes(x = Plants, y = MeansTorPval)) +
  geom_line(colour = "black",size = 1) + ylab("Mean P-value") + 
  xlab("Number of plants") + geom_point(size = 2) + coord_cartesian(ylim = c(0.002, 0.015))+
  scale_y_continuous(breaks = seq(from = 0.002, to = 0.015, by = 0.002))+
  scale_x_reverse(breaks = seq(from = 10, to = 20, by = 2))
plotTorPVal + ng1

# Generate distribution of slopes
S1T <- replicate(1000, FunSlope(KT_plants_TorB, size = 20))
S2T <- replicate(1000, FunSlope(KT_plants_TorB, size = 19))
S3T <- replicate(1000, FunSlope(KT_plants_TorB, size = 18))
S4T <- replicate(1000, FunSlope(KT_plants_TorB, size = 17))
S5T <- replicate(1000, FunSlope(KT_plants_TorB, size = 16))
S6T <- replicate(1000, FunSlope(KT_plants_TorB, size = 15))
S7T <- replicate(1000, FunSlope(KT_plants_TorB, size = 14))
S8T <- replicate(1000, FunSlope(KT_plants_TorB, size = 13))
S9T <- replicate(1000, FunSlope(KT_plants_TorB, size = 12))
S10T <- replicate(1000, FunSlope(KT_plants_TorB, size = 11))
S11T <- replicate(1000, FunSlope(KT_plants_TorB, size = 10))

#Bind results into table#
SlopeTor <- cbind(S1T, S2T, S3T, S4T, S5T, S6T, S7T, S8T, S9T, S10T, S11T)
SlopeTor <- as.data.frame(SlopeTor)
Plants <- c(20, 19, 18, 17, 16, 15, 14, 13, 12 ,11 ,10)
MeansTorSlope <- colMeans(SlopeTor)
SdTorSlopePop <- sapply(SlopeTor, function(cl) sds = sd(cl))
TorPlotSlope <- as.data.frame(cbind(SdTorSlopePop, Plants, MeansTorSlope))

#Plot number of plants vs. mean slope#
plotTorSlope<-ggplot(TorPlotSlope, aes(x = Plants, y = MeansTorSlope)) +
  geom_line(colour = "black",size = 1) + ylab("Mean Slope") + 
  xlab("Number of plants") + geom_point(size = 2) + 
  coord_cartesian(ylim=c(0.002, 0.015)) +
  scale_y_continuous(breaks = seq(from = 0.002, to = 0.015, by = 0.002))+
  scale_x_reverse(breaks = seq(from = 10, to = 20, by = 2))
plotTorSlope + ng1

#Plot number of plants vs. StDev in slope#
plotTorSlopeSD<-ggplot(TorPlotSlope, aes(x = Plants, y = SdTorSlopePop)) +
  geom_line(colour = "black",size = 1) + ylab("Mean Slope") + 
  xlab("Number of plants") + geom_point(size = 2) +
  # coord_cartesian(ylim=c(0.002, 0.015))+
  # scale_y_continuous(breaks = seq(from = 0.002, to = 0.015, by = 0.002))+
  scale_x_reverse(breaks = seq(from = 10, to = 20, by = 2))
plotTorSlopeSD + ng1

### SAMPLING POPULATIONS ###

# Pvals
P1T <- replicate(1000, FunPvalPop(KT_pops_TorB, size = 50))    #50 populations, with replacement
P2T <- replicate(1000, FunPvalPop(KT_pops_TorB, size = 45))    #45 populations, with replacement  #80% of plants per population, without replacement
P3T <- replicate(1000, FunPvalPop(KT_pops_TorB, size = 40))    #40 populations, with replacement
P4T <- replicate(1000, FunPvalPop(KT_pops_TorB, size = 35))    #Etc...
P5T <- replicate(1000, FunPvalPop(KT_pops_TorB, size = 30))
P6T <- replicate(1000, FunPvalPop(KT_pops_TorB, size = 25))
P7T <- replicate(1000, FunPvalPop(KT_pops_TorB, size = 20))
P8T <- replicate(1000, FunPvalPop(KT_pops_TorB, size = 15))
P9T <- replicate(1000, FunPvalPop(KT_pops_TorB, size = 10))

#Bind results into table#
PvalTorPop <- cbind(P1T, P2T, P3T, P4T, P5T, P6T, P7T, P8T, P9T)
PvalTorPop <- as.data.frame(PvalTorPop)
Pops <- c(50, 45, 40, 35, 30, 25, 20, 15, 10)
MeansTorPvalPop <- colMeans(PvalTorPop)
TorPlotPvalPop <- as.data.frame(cbind(Pops, MeansTorPvalPop))

#Plot Means P-value vs. number of populations for Toronto#
plotTorPValPop <- ggplot(TorPlotPvalPop, aes(x = Pops, y = MeansTorPvalPop)) +
  geom_line(colour = "black",size = 1) + ylab("Mean P-value") + 
  xlab("Number of populations") + geom_point(size = 2) + 
  coord_cartesian(ylim = c(0, 0.102)) +
  scale_y_continuous(breaks=seq(from = 0, to = 0.1, by = 0.02)) +
  scale_x_reverse(breaks = seq(from = 10, to = 50, by = 5))
plotTorPValPop + ng1

# Slopes
S1T <- replicate(1000, FunSlopePop(KT_pops_TorB, size = 50))
S2T <- replicate(1000, FunSlopePop(KT_pops_TorB, size = 45))
S3T <- replicate(1000, FunSlopePop(KT_pops_TorB, size = 40))
S4T <- replicate(1000, FunSlopePop(KT_pops_TorB, size = 35))
S5T <- replicate(1000, FunSlopePop(KT_pops_TorB, size = 30))
S6T <- replicate(1000, FunSlopePop(KT_pops_TorB, size = 25))
S7T <- replicate(1000, FunSlopePop(KT_pops_TorB, size = 20))
S8T <- replicate(1000, FunSlopePop(KT_pops_TorB, size = 15))
S9T <- replicate(1000, FunSlopePop(KT_pops_TorB, size = 10))

#Bind results into table#
SlopeTorPop <- cbind(S1T, S2T, S3T, S4T, S5T, S6T, S7T, S8T, S9T)
SlopeTorPop <- as.data.frame(SlopeTorPop)
Pops <- c(50, 45, 40, 35, 30, 25, 20, 15, 10)
MeansTorSlopePop <- colMeans(SlopeTorPop)
SdTorSlopePop <- sapply(SlopeTorPop, function(cl) sds = sd(cl))
TorPlotSlopePop <- as.data.frame(cbind(SdTorSlopePop, Pops, MeansTorSlopePop))

#Plot mean slope vs. number of populations#
plotTorSlopePop <- ggplot(TorPlotSlopePop, aes(x = Pops, y = MeansTorSlopePop)) +
  geom_line(colour = "black", size = 1) + ylab("Mean Slope") + 
  xlab("Number of plants") + geom_point(size = 2) + 
  coord_cartesian(ylim = c(0.002,0.015)) +
  scale_y_continuous(breaks=seq(from = 0.002, to = 0.015,by = 0.002))+
  scale_x_reverse(breaks = seq(from = 10,to = 50, by = 5))
plotTorSlopePop + ng1

#Plot Sd slope vs. number of populations#
plotTorSdSlopePop <- ggplot(TorPlotSlopePop, aes(x = Pops, y = SdTorSlopePop)) +
  geom_line(colour = "black",size = 1) + ylab("StDev Slope") + 
  xlab("Number of plants") + geom_point(size = 2) + 
  coord_cartesian(ylim = c(0.0010,0.0035)) +
  scale_y_continuous(breaks = seq(from = 0.001, to = 0.0035, by = 0.00035))+
  scale_x_reverse(breaks = seq(from = 10, to = 50, by = 5))
plotTorSdSlopePop + ng1
