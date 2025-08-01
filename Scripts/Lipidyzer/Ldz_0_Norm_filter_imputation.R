##Lipidyzer data on primary NSC culture
rm(list=ls())
library(tidyverse)

setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Pre-processing_functions.R")

#### Append sample name to data frame ####===============================================================
# Upload data
all <- read.csv("./Input_Data/a_L_species_Conc_BlkNorm.csv",
                stringsAsFactors = F) 
#exclude QC + Blk, import DNA normalization factor
all.df <- all %>%
  select(., -matches("QC|B"))
#import sample name spreadsheet
key <- read.csv("./Input_Data/sample_key.csv", stringsAsFactors = F)
key <- key %>%
  mutate(., ID = paste0("X", Key))
#Append sample name to data matrix
all_features <- all.df %>%
  rename_at(vars(key$ID), ~key$Samples) %>%
  column_to_rownames(., var = "features")

save(all_features, file = paste0("./Output_Data/Lipidyzer_all_feature.Rdata"))

#### Median concentration normalization ####=============================================================
load("./Output_Data/Lipidyzer_all_feature.Rdata")
Med_ind <- apply(all_features, 2, median, na.rm=T)
Med_norm_factor <- as.numeric(scale(Med_ind, center = F, scale = median(Med_ind)))

Med_norm.ldz <- sweep(all_features, 2, Med_norm_factor, "/")
N_med <- apply(Med_norm.ldz, 2, median, na.rm=T)

save(Med_norm.ldz, file = paste0("./Output_Data/Ldz_Med_conc_normed_all_features.Rdata"))

#### Filter and Imputation####===========================================================================
##filter lipid should be detected at least in half of the samples (NA <10)
load("./Output_Data/Ldz_Med_conc_normed_all_features.Rdata")
Norm_filter <- filter.Ldz(Med_norm.ldz, 21) ##511 lipids meet this criteria

##Imputation to replace missing value
Norm_filter.ipt <- impt.Ldz(Norm_filter)

##Check the number of values that were replace by imputation = number of missing value to begin with
impt.check.Ldz(Norm_filter, Norm_filter.ipt)

Ldz.norm.impt <- Norm_filter.ipt
save(Ldz.norm.impt, file = paste0("./Output_Data/Ldz_511_lipids_log2_MedConc_norm_imputed.Rdata"))

##Also create a dataframe with values transform back to raw concentration####
ldz.raw_conc <- 2^(Ldz.norm.impt)
save(ldz.raw_conc, file = paste0("./Output_data/Ldz_511_lipids_backtoraw.conc_MedConc_norm_imputed.Rdata"))


