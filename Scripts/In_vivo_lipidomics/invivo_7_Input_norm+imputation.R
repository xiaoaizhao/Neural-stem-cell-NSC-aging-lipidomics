
##Steps:
##1. Apply median concentration normalization to adjust for variation in input material. Use concentration from 121 lipids with concentration, 
##and then apply the same normalization factor to normalize 121 lipids with concentration, as well as all 130 lipid (including 9 without concentration)
##2. Further filter to only keep lipids that are detected in at least half of the samples (all lipids meet this criteria)
##3. Remove NA value by imputation - replace with imputed value from the bottom 10% intensity of each lipids
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load(file = "./Output_Data/Invivo_cell_121_lipid_conc.Rdata")

conc.only <- lipid.w.conc %>%
  ungroup() %>%
  select(., -c("Class", "Ion")) %>%
  column_to_rownames(., var = "LipidIon.FA")

##get normalization factor
MedConc <- apply(conc.only, 2, median, na.rm=T)
Med_ind <- as.numeric(scale(MedConc, center = F, scale = median(MedConc)))

##normalized 121 lipids with concentration
Med_norm_121 <- sweep(conc.only, 2, Med_ind, "/")
newMed <- apply(Med_norm_121, 2, median, na.rm=T) ##same median concentration after norm

##normalized all 130 lipids (121 lipids with concentration +9 without concentration)
load(file = "./Output_Data/Invivo_Cell_All_121_conc_9_no_conc.Rdata")
conc.w.noConc <- Cell.lipid.all %>%
  ungroup() %>%
  select(., -c("Class", "Ion")) %>%
  column_to_rownames(., var = "LipidIon.FA")
Med_norm_130 <- sweep(conc.w.noConc, 2, Med_ind, "/")
newMed130 <- apply(Med_norm_130, 2, median, na.rm=T) ##median after norm is NOT identical, this makes sense becuase normalization factor was generated using 121 lipids with concentration only

###Imputation to remove 0 value####################################################################################
###lipid with concentration first
##filter lipid should be detected at least in half of the samples (NA <6)
Norm_filter <- filter.func(Med_norm_121, 12) ##121, all lipids meet this criteria

##Imputation to replace missing value
Norm_filter.ipt <- impt.func(Med_norm_121)

##Check the number of values that were replace by imputation = number of missing value to begin with
impt.check(Med_norm_121, Norm_filter.ipt)

Invivo.Impt_norm_conc_all <- Norm_filter.ipt
save(Invivo.Impt_norm_conc_all, file = paste0("./Output_Data/Invivo_Norm_Impt_log2_conc_121_lipids.Rdata"))

Invivo.Impt.norm.conc.raw <- 2^Invivo.Impt_norm_conc_all %>% 
  rownames_to_column(var = "LipidIon") %>% 
  pivot_longer(-LipidIon, values_to = "Concentration", names_to = "Samples") %>% 
  mutate(Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>% 
  mutate(Concentration = ifelse(Class == "LPC", Concentration/50, Concentration)) %>% 
  pivot_wider(names_from = "Samples", values_from = "Concentration")

save(Invivo.Impt.norm.conc.raw, file = paste0("./Output_Data/Invivo_Norm_Impt_RAW.conc_121_lipids.Rdata"))
###Imputation to remove 0 value####################################################################################
###lipid with concentration + without concentration
##filter lipid should be detected at least in half of the samples (NA <6)
Norm_filter.all <- filter.func(Med_norm_130, 12) ##130, all lipids meet this criteria

##Imputation to replace missing value
Norm_filter.ipt.all <- impt.func(Med_norm_130)

##Check the number of values that were replace by imputation = number of missing value to begin with
impt.check(Med_norm_130, Norm_filter.ipt.all)

Invivo.Impt_norm_conc_no_conc_all <- Norm_filter.ipt.all
save(Invivo.Impt_norm_conc_no_conc_all, file = paste0("./Output_Data/Invivo_Norm_Impt_log2_all130_lipids.Rdata"))

##Also create a dataframe with values transform back to raw intensity####
Invivo.raw_int.invivo <- 2^(Invivo.Impt_norm_conc_no_conc_all)
save(Invivo.raw_int.invivo, file = paste0("./Output_data/Invivo_Norm_Impt_backtoraw_all130_lipids.Rdata"))
