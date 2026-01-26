### Normalization and imputation

setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)

source("./Scripts/Function_scripts/Pre-processing_functions.R")
load(file = "./Output_data/Invivo_29.lipids.conc.byMsD.Rdata")

## Normalization and imputation should exclude all QC data, include actual samples only, remove QC samples and prep blank
conc.only <- Invivo.lipid.w.conc.MsD.29peak %>%
  ungroup() %>%
  filter(!grepl("d7|d9", Metabolite))

Invivo.for.norm <- conc.only %>% 
  column_to_rownames(., var = "Metabolite") %>% 
  arrange(row.names(.)) %>% 
  select(order(colnames(.)))

##get normalization factor
MedConc <- apply(Invivo.for.norm, 2, median, na.rm=T)
Med_ind <- as.numeric(scale(MedConc, center = F, scale = median(MedConc)))

##normalized 48 lipids with concentration
Med_norm_29 <- sweep(Invivo.for.norm, 2, Med_ind, "/")
newMed <- apply(Med_norm_29, 2, median, na.rm=T) ##same median concentration after norm
unique(newMed) #0.2369455


## ====Imputation lipid with concentration ====
###Imputation to remove 0 value####################################################################################

##filter lipid should be detected at least in half of the samples (NA <6)
Norm_filter <- filter.func(Med_norm_29, 12) ##29 out of 29 lipids meet this criteria

##Imputation to replace missing value
Norm_filter.ipt <- impt.func(Norm_filter)


##Check the number of values that were replace by imputation = number of missing value to begin with
Norm_filter <- Norm_filter %>% 
  replace(is.na(.), 0)
impt.check(Norm_filter, Norm_filter.ipt)

Invivo.log2.impt.norm.conc.MsD.good.peak <- Norm_filter.ipt
save(Invivo.log2.impt.norm.conc.MsD.good.peak, file = paste0("./Output_data/Invivo_MsD.Norm_Impt_log2_conc_29.lipids.Rdata"))
