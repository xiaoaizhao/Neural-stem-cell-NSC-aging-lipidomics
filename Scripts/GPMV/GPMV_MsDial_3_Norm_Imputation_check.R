### Normalization and imputation


setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)

source("./Scripts/Function_scripts/Pre-processing_functions.R")
load(file = "./Output_data/GPMV_291.quant.lipids.Rdata")

## Normalization and imputation should exclude all QC data, include actual samples only, remove QC samples and prep blank
conc.only <- GPMV.lipid.w.conc %>%
  ungroup() %>%
  filter(!grepl("d7|d9", Metabolite.name))

GPMV.for.norm <- conc.only %>% 
  column_to_rownames(., var = "Metabolite.name") %>% 
  arrange(row.names(.)) %>% 
  select(order(colnames(.)))

##get normalization factor
MedConc <- apply(GPMV.for.norm, 2, median, na.rm=T)
Med_ind <- as.numeric(scale(MedConc, center = F, scale = median(MedConc)))

##normalized lipids with concentration
Med_norm_291 <- sweep(GPMV.for.norm, 2, Med_ind, "/")
newMed <- apply(Med_norm_291, 2, median, na.rm=T) ##same median concentration after norm
unique(newMed) #0.2359657


## ====Imputation lipid with concentration ====
###Imputation to remove 0 value####################################################################################

##filter lipid should be detected at least in half of the samples (NA <6)
Norm_filter <- filter.func(Med_norm_291, 16) ##291 out of 291 lipids meet this criteria

##Imputation to replace missing value
Norm_filter.ipt <- impt.func(Norm_filter)


##Check the number of values that were replace by imputation = number of missing value to begin with
Norm_filter <- Norm_filter %>% 
  replace(is.na(.), 0)
impt.check(Norm_filter, Norm_filter.ipt)

GPMV.log2.impt.norm.conc.MsD <- Norm_filter.ipt
save(GPMV.log2.impt.norm.conc.MsD, file = paste0("./Output_data/GPMV_MsD.Norm_Impt_log2_conc_291_lipids.Rdata"))


GPMV.impt.norm.raw.conc.lpd.MsD <- 2^GPMV.log2.impt.norm.conc.MsD
save(GPMV.impt.norm.raw.conc.lpd.MsD, file = "./Output_Data/GPMV_MsD.Raw.conc.Norm_Impt_291_lipids.Rdata")
