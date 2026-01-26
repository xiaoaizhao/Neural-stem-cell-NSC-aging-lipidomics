### Normalization and imputation

setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)

source("./Scripts/Function_scripts/Pre-processing_functions.R")
load(file = "./Output_data/Exp2_278.quant.lipids.byMsD.Rdata")

## Normalization and imputation should exclude all QC data, include actual samples only, remove QC samples and prep blank
conc.only <- Exp2.lipid.w.conc.MsD %>%
  ungroup() %>%
  filter(!grepl("d7|d9", Metabolite.name))

Exp2.for.norm <- Exp2.lipid.w.conc.MsD %>% 
  rowwise() %>%
  mutate(LipidIon = Metabolite.name) %>%
  relocate(LipidIon, .after = Metabolite) %>%
  column_to_rownames(., var = "LipidIon") %>%
  select(-c(Metabolite.name, Metabolite))

##get normalization factor
MedConc <- apply(Exp2.for.norm, 2, median, na.rm=T)
Med_ind <- as.numeric(scale(MedConc, center = F, scale = median(MedConc)))

##normalized 278 lipids with concentration
Med_norm_278 <- sweep(Exp2.for.norm, 2, Med_ind, "/")
newMed <- apply(Med_norm_278, 2, median, na.rm=T) ##same median concentration after norm
unique(newMed) #0.09998704


## ====Imputation lipid with concentration ====
###Imputation to remove 0 value####################################################################################

##filter lipid should be detected at least in half of the samples (NA <6)
Norm_filter <- filter.func(Med_norm_278, 48) ##278 out of 287 lipids meet this criteria

##Imputation to replace missing value
Norm_filter.ipt <- impt.func(Norm_filter)

##Check the number of values that were replace by imputation = number of missing value to begin with
Norm_filter <- Norm_filter %>% 
  replace(is.na(.), 0)
impt.check(Norm_filter, Norm_filter.ipt)

Exp2.log2.impt.norm.conc.MsD <- Norm_filter.ipt
save(Exp2.log2.impt.norm.conc.MsD, file = paste0("./Output_data/Exp2_MsD.Norm_Impt_log2_conc_278_lipids.Rdata"))

