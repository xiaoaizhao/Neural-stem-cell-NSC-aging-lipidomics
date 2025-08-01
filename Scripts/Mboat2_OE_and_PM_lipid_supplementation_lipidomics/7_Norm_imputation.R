### Processing script 7, _Mboat2_ overexpression and young plasma membrane lipid supplementation lipidomics 
#1. input norm + imputation
# 275 lipids with conc, 48 lipids without conc
setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(reshape2)
library(stringi)

source("./Scripts/Function_scripts/Pre-processing_functions.R")
load(file = "./Output_data/M2PM_275_conc_lipids.Rdata")

## Normalization and imputation should exclude all QC data, include actual samples only
conc.only <- M2PM.lipid.w.conc %>%
  ungroup() %>%
  select(., -contains("QC")) %>%
  column_to_rownames(., var = "LipidIon")

##get normalization factor
MedConc <- apply(conc.only, 2, median, na.rm=T)
Med_ind <- as.numeric(scale(MedConc, center = F, scale = median(MedConc)))

##normalized 275 lipids with concentration
Med_norm_275 <- sweep(conc.only, 2, Med_ind, "/")
newMed <- apply(Med_norm_275, 2, median, na.rm=T) ##same median concentration after norm
unique(newMed)

##normalized all 323 lipids (275 lipids with concentration +48 without concentration)
load(file = "./Output_data/M2PM_275_conc_48_noconc_lipids.Rdata")

conc.w.noConc <- M2PM_lpd.all %>% #323 lipids total
  ungroup() %>%
  select(., -contains("QC")) %>%
  column_to_rownames(., var = "LipidIon")
Med_norm_323 <- sweep(conc.w.noConc, 2, Med_ind, "/")
newMed323 <- apply(Med_norm_323, 2, median, na.rm=T) ##median after norm is NOT identical, this makes sense becuase normalization factor was generated using 447 lipids with concentration only
unique(newMed323)
###Imputation to remove 0 value####################################################################################
###lipid with concentration first
##filter lipid should be detected at least in half of the samples (NA <6)
Norm_filter <- filter.func(Med_norm_275, 40) ##273 out of 275 lipids meet this criteria, all pass

##Imputation to replace missing value
Norm_filter.ipt <- impt.func(Norm_filter)

##Check the number of values that were replace by imputation = number of missing value to begin with
Norm_filter <- Norm_filter %>% 
  replace(is.na(.), 0)
impt.check(Norm_filter, Norm_filter.ipt)

M2PM.Impt_norm_conc_all <- Norm_filter.ipt
save(M2PM.Impt_norm_conc_all, file = paste0("./Output_data/M2PM_Norm_Impt_log2_conc_275_lipids.Rdata"))

###Imputation to remove 0 value####################################################################################
###lipid with concentration + without concentration
##filter lipid should be detected at least in half of the samples 
Norm_filter.all <- filter.func(Med_norm_323, 40) ##321 out of 323 lipids meet this criteria, all pass

##Imputation to replace missing value
Norm_filter.ipt.all <- impt.func(Norm_filter.all)

##Check the number of values that were replace by imputation = number of missing value to begin with
Norm_filter.all <- Norm_filter.all %>% 
  replace(is.na(.), 0)
impt.check(Norm_filter.all, Norm_filter.ipt.all)

M2PM.Impt_norm_conc_no_conc_all <- Norm_filter.ipt.all
save(M2PM.Impt_norm_conc_no_conc_all, file = paste0("./Output_data/M2PM_Norm_Impt_log2_conc_275_noconc_48_lipids.Rdata"))

##Also create a dataframe with values transform back to raw intensity####
M2PM.raw.all<- 2^(M2PM.Impt_norm_conc_no_conc_all)
save(M2PM.raw.all, file = paste0("./Output_data/M2PM_backtoRAW_Norm_Impt_log2_conc_275_noconc_48_lipids.Rdata"))
