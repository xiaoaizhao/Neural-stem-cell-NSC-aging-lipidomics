### Normalization and imputation

setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)

source("./Scripts/Function_scripts/Pre-processing_functions.R")
load(file = "./Output_data/M2PM_397.quant.lipids.byMsD.Rdata")

## Normalization and imputation should exclude all QC data, include actual samples only, remove QC samples and prep blank
conc.only <- M2PM.lipid.w.conc.MsD %>%
  ungroup() %>%
  filter(!grepl("d7|d9", Metabolite.name))

M2PM.for.norm <- conc.only %>% 
  column_to_rownames(., var = "Metabolite.name") %>% 
  arrange(row.names(.)) %>% 
  select(order(colnames(.)))

##get normalization factor
MedConc <- apply(M2PM.for.norm, 2, median, na.rm=T)
Med_ind <- as.numeric(scale(MedConc, center = F, scale = median(MedConc)))

##normalized 48 lipids with concentration
Med_norm_397 <- sweep(M2PM.for.norm, 2, Med_ind, "/")
newMed <- apply(Med_norm_397, 2, median, na.rm=T) ##same median concentration after norm
unique(newMed) #0.06700575


## ====Imputation lipid with concentration ====
###Imputation to remove 0 value####################################################################################

##filter lipid should be detected at least in half of the samples (NA <6)
Norm_filter <- filter.func(Med_norm_397, 40) ##397 out of 397 lipids meet this criteria

##Imputation to replace missing value
Norm_filter.ipt <- impt.func(Norm_filter)


##Check the number of values that were replace by imputation = number of missing value to begin with
Norm_filter <- Norm_filter %>% 
  replace(is.na(.), 0)
impt.check(Norm_filter, Norm_filter.ipt)

M2PM.log2.impt.norm.conc.MsD <- Norm_filter.ipt
save(M2PM.log2.impt.norm.conc.MsD, file = paste0("./Output_data/M2PM_MsD.Norm_Impt_log2_conc_397_lipids.Rdata"))


M2PM.MsD.lpd.conc.all<- 2^(M2PM.log2.impt.norm.conc.MsD)
save(M2PM.MsD.lpd.conc.all, file = paste0("./Output_data/M2PM_backtoRAW_MsD_Norm_Impt_397_lipids.Rdata"))
