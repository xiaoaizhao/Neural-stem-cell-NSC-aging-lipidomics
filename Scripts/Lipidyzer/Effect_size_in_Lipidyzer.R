##Effect size calculation for Lipidzyer data
##Use Hedge's g effect size calculation
##calculated effect size both for act vs. qui cells and old vs. young in quiescent cell only
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")
##Effect size calculation for lipidyzer data####
load(file = "./Output_Data/Ldz_511_lipids_backtoraw.conc_MedConc_norm_imputed.Rdata")

eff.size.dyzer <- ldz.raw_conc %>%
  rownames_to_column(., var = "Lipid") %>%
  pivot_longer(-Lipid, names_to = "Samples", values_to = "Concentration") %>%
  mutate(., Age = ifelse(grepl("Y", Samples, ignore.case = F), "Young", "Old")) %>%
  mutate(., CellType = ifelse(grepl("act", Samples), "Activated", "Quiescent"))

##Effect size between quiescent and activated cells####
Lipidyzer.cell.es.g <- es.g.func.celltype(eff.size.dyzer, Lipid, CellType, Concentration, Samples) 

save(Lipidyzer.cell.es.g, file = paste0("./Output_Data/Ldz_Ef_Size_celltype.Rdata"))

##Effect size between old and young quiescent cells####
Qui.dyzer <- eff.size.dyzer %>%
  filter(., CellType == "Quiescent")
Lipidyzer.Age.es.g <- es.g.func(Qui.dyzer, Lipid, Age, Concentration, Samples) 

save(Lipidyzer.Age.es.g, file = paste0("./Output_Data/Ldz_Ef_Size_Age_qNSC.Rdata"))
