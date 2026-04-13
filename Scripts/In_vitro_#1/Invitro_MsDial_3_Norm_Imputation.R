### Normalization and imputation

setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)

source("./Scripts/Function_scripts/Pre-processing_functions.R")
load(file = "./Output_data/Invitro_404.quant.lipids.Smpl.Name.O8aNSC.adj.Rdata")

## Normalization and imputation should exclude all QC data, include actual samples only, remove QC samples and prep blank
conc.only <- Invitro.lpd.conc.O8aNSC.adj %>%
  ungroup() %>%
  filter(!grepl("d7|d9", Metabolite.name)) %>% 
  select(-Metabolite)

df.for.norm <- conc.only %>% 
  column_to_rownames(., var = "Metabolite.name") %>% 
  arrange(row.names(.)) %>% 
  select(order(colnames(.)))


##get normalization factor
MedConc <- apply(df.for.norm, 2, median, na.rm=T)
Med_ind <- as.numeric(scale(MedConc, center = F, scale = median(MedConc)))

##normalized 404 lipids with concentration
Med_norm_404 <- sweep(df.for.norm, 2, Med_ind, "/")
newMed <- apply(Med_norm_404, 2, median, na.rm=T) ##same median concentration after norm
unique(newMed) #0.07669466


## ====Imputation lipid with concentration ====
###Imputation to remove 0 value####################################################################################

##filter lipid should be detected at least in half of the samples
Norm_filter <- filter.func(Med_norm_404, 20) ##404 out of 404 lipids meet this criteria

##Imputation to replace missing value
Norm_filter.ipt <- impt.func(Norm_filter)


##Check the number of values that were replace by imputation = number of missing value to begin with
Norm_filter <- Norm_filter %>% 
  replace(is.na(.), 0)
impt.check(Norm_filter, Norm_filter.ipt)

Invitro.log2.impt.norm.conc.MsD <- Norm_filter.ipt
save(Invitro.log2.impt.norm.conc.MsD, file = paste0("./Output_data/Invitro_MsD.Norm_Impt_log2_conc_404_lipids.Rdata"))
