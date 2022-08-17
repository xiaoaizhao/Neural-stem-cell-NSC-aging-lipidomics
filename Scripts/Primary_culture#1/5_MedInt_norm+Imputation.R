##Objective: Median intensity normalization + filter lipid + imputation to remove missing value

##Steps:
##1. Apply median intensity normalization to adjust for variation in input material
##2. Further filter to only keep lipids that are detected in at least half of the samples (all lipids meet this criteria)
##3. Remove NA value by imputation - replace with imputated value from the bottom 10% intensity of each lipids
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load(file = "./Output_data/Spike-in_norm_373_lipid_w_11Qui_avg.Rdata")

##normalized by median intensity of each sample
dftonorm <- tg_norm_lipid_373_w_11_avg %>%
  column_to_rownames(., var = "LipidIon")

MedConc <- apply(dftonorm, 2, median, na.rm=T)
Med_ind <- as.numeric(scale(MedConc, center = F, scale = median(MedConc)))

##normalized 373 lipids with intensity
Med_norm_373 <- sweep(dftonorm, 2, Med_ind, "/")
newMed <- apply(Med_norm_373, 2, median, na.rm=T) 
unique(newMed) ##[1] 393520.9 - same median intensity after norm

###Imputation to remove 0 value####################################################################################
##filter lipid should be detected at least in half of the samples (NA <12)
Norm_filter <- filter.func(Med_norm_373, 24) ##373, all lipids meet this criteria

##Imputation to replace missing value
Norm_filter.ipt <- impt.func(Med_norm_373)

##Check the number of values that were replace by imputation = number of missing value to begin with
impt.check(Med_norm_373, Norm_filter.ipt)

Impt_norm_373all <- Norm_filter.ipt
save(Impt_norm_373all, file = paste0("./Output_Data/Spike-in_norm_MedNorm_all_373_lipid.Rdata"))


##Also create a data frame with values transform back to raw intensity####
raw_int <- 2^(Impt_norm_373all)
save(raw_int, file = paste0("./Output_data/Spike-in_norm_Mednorm_all_373_lipids_back_to_raw_int.Rdata"))
