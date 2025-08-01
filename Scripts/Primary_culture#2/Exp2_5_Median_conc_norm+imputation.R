
##Steps:
##1. Apply median concentration normalization to adjust for variation in input material. Use concentration from 613 lipids with concentration, 
##and then apply the same normalization factor to normalize 613 lipids with concentration, as well as all 694 lipid (including 81 without concentration)
##2. Further filter to only keep lipids that are detected in at least half of the samples (all lipids meet this criteria)
##3. Remove NA value by imputation - replace with imputed value from the bottom 10% intensity of each lipids
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())


####Apply median concentration normalization to adjust for variation in input material##################################################
load(file = "./Output_Data/Exp2_613_lipids_w_conc.Rdata") ##613 lipids with concentration

conc.only <- exp2.lipid.w.conc %>%
  ungroup(Class, Ion) %>%
  select(., matches("Lipid|Y|O", ignore.case = F)) %>%
  column_to_rownames(., var = "LipidIon")

MedConc <- apply(conc.only, 2, median, na.rm=T)
Med_ind <- as.numeric(scale(MedConc, center = F, scale = median(MedConc)))

#### Get norm factor for later FFA analysis
norm_mtx <- as.data.frame(MedConc)
norm_mtx$ind <-  as.numeric(scale(norm_mtx$MedConc, center = F, scale = median(MedConc)))
save(norm_mtx, file = paste0("./Output_Data/Exp2_Norm_factor_from_613lpds.Rdata"))

##normalized 613 lipids with concentration
Med_norm_613 <- sweep(conc.only, 2, Med_ind, "/")
newMed <- apply(Med_norm_613, 2, median, na.rm=T) ##same median concentration after norm 0.0438
save(Med_norm_613, file = paste0("./Output_Data/Exp2_Med_conc_norm_613_lipids.Rdata"))

##normalized all 694 lipids (613 lipids with concentration +81 without concentration)
load(file = "./Output_Data/Exp2_All_613_conc_81_no_conc.Rdata")
all.694 <- exp2.lipid.all %>%
  ungroup(Class, Ion) %>%
  select(., matches("Lipid|Y|O", ignore.case = F)) %>%
  column_to_rownames(., var = "LipidIon")

Med_norm_694 <- sweep(all.694, 2, Med_ind, "/")
newMed694 <- apply(Med_norm_694, 2, median, na.rm=T) 
save(Med_norm_694, file = paste0("./Output_Data/Exp2_Med_conc_norm_694_lipids_(613conc_81int).Rdata"))

###Imputation to remove 0 value####################################################################################
###lipid with concentration first
##filter lipid should be detected at least in half of the samples (NA <24)
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Exp2_Med_conc_norm_613_lipids.Rdata")
Norm_filter <- filter.func(Med_norm_613, 48) ##612, 1 lipid does not meet the criteria

##Imputation to replace missing value
Norm_filter.ipt <- impt.Ldz(Norm_filter)

##Check the number of values that were replace by imputation = number of missing value to begin with
impt.check.Ldz(Norm_filter, Norm_filter.ipt)

Exp2_Impt_norm_conc_all <- Norm_filter.ipt
save(Exp2_Impt_norm_conc_all, file = paste0("./Output_Data/Exp2_Norm_Impt_log2_conc_612_lipids.Rdata"))

##Also create a dataframe with values transform back to raw intensity####
raw_conc.exp2 <- 2^(Exp2_Impt_norm_conc_all)
save(raw_conc.exp2, file = paste0("./Output_data/Exp2_Norm_Impt_backtoraw_conc612_lipids.Rdata"))

###Imputation to remove 0 value####################################################################################
###lipid with concentration + without concentration
##filter lipid should be detected at least in half of the samples (NA <24)
load("./Output_Data/Exp2_Med_conc_norm_694_lipids_(613conc_81int).Rdata")
Norm_filter.all <- filter.func(Med_norm_694, 48) ##693, 1 lipid does not meet the criteria

##Imputation to replace missing value
Norm_filter.ipt.all <- impt.Ldz(Norm_filter.all)

##Check the number of values that were replace by imputation = number of missing value to begin with
impt.check.Ldz(Norm_filter.all, Norm_filter.ipt.all)

Exp2_Impt_norm_conc_no_conc_all <- Norm_filter.ipt.all
save(Exp2_Impt_norm_conc_no_conc_all, file = paste0("./Output_Data/Exp2_Norm_Impt_log2_all693_lipids.Rdata"))

##Also create a dataframe with values transform back to raw intensity####
raw_int.exp2 <- 2^(Exp2_Impt_norm_conc_no_conc_all)
save(raw_int.exp2, file = paste0("./Output_data/Exp2_Norm_Impt_backtoraw_all693_lipids.Rdata"))
