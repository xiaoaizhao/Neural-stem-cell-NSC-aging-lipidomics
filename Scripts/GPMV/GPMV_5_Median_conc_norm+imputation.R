##Steps:
##1. Apply median concentration normalization to adjust for variation in input material. Use concentration from 483 lipids with concentration, 
##and then apply the same normalization factor to normalize 483 lipids with concentration, as well as all 565 lipid (including 82 without concentration)
##2. Further filter to only keep lipids that are detected in at least half of the samples (all lipids meet this criteria)
##3. Remove NA value by imputation - replace with imputed value from the bottom 10% intensity of each lipids
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

####Apply median concentration normalization to adjust for variation in input material##################################################
load(file = "./Output_Data/GPMV_483lipid_concentration.Rdata") ##483 lipids with concentration

conc.only <- GPMV.w.conc %>%
  ungroup(Class, Ion) %>%
  select(., matches("Lipid|Y|O", ignore.case = F)) %>%
  column_to_rownames(., var = "LipidIon")

MedConc <- apply(conc.only, 2, median, na.rm=T)
Med_ind <- as.numeric(scale(MedConc, center = F, scale = median(MedConc)))

##normalized 483 lipids with concentration
Med_norm_483 <- sweep(conc.only, 2, Med_ind, "/")
newMed <- apply(Med_norm_483, 2, median, na.rm=T) ##same median concentration after norm 0.159
save(Med_norm_483, file = paste0("./Output_Data/GPMV_Med_conc_norm_483_lipids.Rdata"))

##normalized all 565 lipids (483 lipids with concentration +82 without concentration)
load(file = "./Output_Data/GPMV_All_483_conc_82_no_conc.Rdata")
all.565 <- GPMV.565.lpd.all %>%
  ungroup(Class, Ion) %>%
  select(., matches("Lipid|Y|O", ignore.case = F)) %>%
  column_to_rownames(., var = "LipidIon")

Med_norm_565 <- sweep(all.565, 2, Med_ind, "/")
newMed565 <- apply(Med_norm_565, 2, median, na.rm=T) 
save(Med_norm_565, file = paste0("./Output_Data/GPMV_Med_conc_norm_565_lipids_(483conc_82int).Rdata"))

###Imputation to remove 0 value####################################################################################
###lipid with concentration first
##filter lipid should be detected at least in half of the samples (NA <8)
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/GPMV_Med_conc_norm_483_lipids.Rdata")
Norm_filter <- filter.func(Med_norm_483, 16) ##all lipid pass

##Imputation to replace missing value
Norm_filter.ipt <- impt.func(Norm_filter)

##Check the number of values that were replace by imputation = number of missing value to begin with
impt.check(Norm_filter, Norm_filter.ipt)

GPMV_Impt_norm_conc_all <- Norm_filter.ipt
save(GPMV_Impt_norm_conc_all, file = paste0("./Output_Data/GPMV_Norm_Impt_log2_conc_483_lipids.Rdata"))

###Imputation to remove 0 value####################################################################################
###lipid with concentration + without concentration
##filter lipid should be detected at least in half of the samples (NA <8)
load("./Output_Data/GPMV_Med_conc_norm_565_lipids_(483conc_82int).Rdata")
Norm_filter.all <- filter.func(Med_norm_565, 16) ##all lipid pass

##Imputation to replace missing value
Norm_filter.ipt.all <- impt.func(Norm_filter.all)

##Check the number of values that were replace by imputation = number of missing value to begin with
impt.check(Norm_filter.all, Norm_filter.ipt.all)

GPMV_Impt_norm_conc_no_conc_all <- Norm_filter.ipt.all
save(GPMV_Impt_norm_conc_no_conc_all, file = paste0("./Output_Data/GPMV_Norm_Impt_log2_all565_lipids.Rdata"))

##Also create a dataframe with values transform back to raw intensity####
raw_int.GPMV <- 2^(GPMV_Impt_norm_conc_no_conc_all)
save(raw_int.GPMV, file = paste0("./Output_data/GPMV_Norm_Impt_backtoraw_all565_lipids.Rdata"))


