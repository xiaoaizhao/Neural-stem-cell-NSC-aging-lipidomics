### Processing script 7, Primary culture #3
#1. input norm + imputation

setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(reshape2)
library(stringi)

source("./Scripts/Function_scripts/Pre-processing_functions.R")
load(file = "./Output_data/Exp3_354_conc_lipidsw2endo.chol_356_O8A-adj.Rdata")

## Normalization and imputation should exclude all QC data, include actual samples only, remove QC samples and prep blank
conc.only <- lipid.w.conc.O8A.adj %>%
  ungroup() %>%
  select(., -contains("QC")) %>%
  select(., -contains("PBS")) %>%
  column_to_rownames(., var = "LipidIon")

##get normalization factor
MedConc <- apply(conc.only, 2, median, na.rm=T)
Med_ind <- as.numeric(scale(MedConc, center = F, scale = median(MedConc)))

##normalized 356 lipids with concentration
Med_norm_356 <- sweep(conc.only, 2, Med_ind, "/")
newMed <- apply(Med_norm_356, 2, median, na.rm=T) ##same median concentration after norm
unique(newMed)

##normalized all 416 lipids (356 (354+2) lipids with concentration +60 without concentration)
load(file = "./Output_data/Exp3_356_conc_60_noconc_lipids_O8A-adj.Rdata")

conc.w.noConc <- Exp3_lpd.all.O8A.adj %>% #416 lipids total
  ungroup() %>%
  select(., -contains("QC")) %>%
  select(., -contains("PBS")) %>%
  column_to_rownames(., var = "LipidIon")
Med_norm_416 <- sweep(conc.w.noConc, 2, Med_ind, "/")
newMed416 <- apply(Med_norm_416, 2, median, na.rm=T) ##median after norm is NOT identical, this makes sense because normalization factor was generated using 351 lipids with concentration only
unique(newMed416)

## ====lipid with concentration first====
###Imputation to remove 0 value####################################################################################

##filter lipid should be detected at least in half of the samples (NA <6)
Norm_filter <- filter.func(Med_norm_356, 20) ##356 out of 356 lipids meet this criteria, all pass

##Imputation to replace missing value
Norm_filter.ipt <- impt.func(Norm_filter)

##Check the number of values that were replace by imputation = number of missing value to begin with
Norm_filter <- Norm_filter %>% 
  replace(is.na(.), 0)
impt.check(Norm_filter, Norm_filter.ipt)

Exp3.Impt_norm_conc_all <- Norm_filter.ipt
save(Exp3.Impt_norm_conc_all, file = paste0("./Output_data/Exp3_Norm_Impt_log2_conc_356_lipids.Rdata"))

## ====lipid with concentration + without concentration====
###Imputation to remove 0 value####################################################################################
##filter lipid should be detected at least in half of the samples 
Norm_filter.all <- filter.func(Med_norm_416, 20) ##416 out of 416 lipids meet this criteria, all pass

##Imputation to replace missing value
Norm_filter.ipt.all <- impt.func(Norm_filter.all)

##Check the number of values that were replace by imputation = number of missing value to begin with
Norm_filter.all <- Norm_filter.all %>% 
  replace(is.na(.), 0)
impt.check(Norm_filter.all, Norm_filter.ipt.all)

Exp3.Impt_norm_conc_no_conc_all <- Norm_filter.ipt.all
save(Exp3.Impt_norm_conc_no_conc_all, file = paste0("./Output_data/Exp3_Norm_Impt_log2_conc_356_noconc_60_lipids.Rdata"))

##Also create a dataframe with values transform back to raw intensity####
raw_int.Exp3<- 2^(Exp3.Impt_norm_conc_no_conc_all)
save(raw_int.Exp3, file = paste0("./Output_data/Exp3_backtoRAW_Norm_Impt_log2_conc_356_noconc_60_lipids.Rdata"))
