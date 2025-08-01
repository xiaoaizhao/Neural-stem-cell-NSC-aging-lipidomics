## _Mboat2_ overexpression and young plasma membrane lipid supplementation lipidomics 
##This is the first step in double bond composition analysis
##Steps:
##1. Fix cholesterol annotation of endogenous cholesterol (one from LipidSearch annotation)
##1. Subset by presence or absence of each given double bond, regardless of how many times (side chains) it was identified
##2. Tally total intensity with every given double bond number in each class.

rm(list = ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load(file = "./Output_data/M2PM_backtoRAW_Norm_Impt_log2_conc_275_noconc_48_lipids.Rdata")

WC <- M2PM.raw.all %>%
  rownames_to_column(var = "LipidIon") %>% 
  mutate(., LipidIon = ifelse(
    grepl("Ch\\+H\\-H2O|Cholesterol", LipidIon),
   "Cholesterol(0:0)+H-H20",
    LipidIon
  ))

WC<- WC %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) #7904/416=19, include both quiescent and activated samples

##Function to get aggregated intensity for each Class:DB combination####
DB_agg <- db.tally(WC, Conc_Int, Sample)
M2PM.DB_by_class <- dplyr::bind_rows(DB_agg)

save(M2PM.DB_by_class, file = paste0("./Output_Data/M2PM_all_sample_DB.Rdata"))


