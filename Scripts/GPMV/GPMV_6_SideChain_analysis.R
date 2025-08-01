##This is the first step in double bond composition analysis
##Steps:
##1. Subset by presence or absence of each given double bond, regardless of how many times (side chains) it was identified
##2. Tally total intensity with every given double bond number in each class.
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Pre-processing_functions.R")

load(file = "./Output_data/GPMV_Norm_Impt_backtoraw_all565_lipids.Rdata")
GPMV.df <- raw_int.GPMV %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) #282

##Function to get aggregated intensity for each Class:DB combination####
DB_agg <- db.tally(GPMV.df, Conc_Int, Sample)
DB_by_class_GPMV <- dplyr::bind_rows(DB_agg) #565

save(DB_by_class_GPMV, file = paste0("./Output_Data/GPMV_DB_by_Class.Rdata"))


#### Side chain analysis on conc lipids only ####
rm(list = ls())
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/GPMV_Norm_Impt_log2_conc_483_lipids.Rdata")

GPMV.df <- 2^GPMV_Impt_norm_conc_all %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) #282

##Function to get aggregated intensity for each Class:DB combination####
DB_agg <- db.tally(GPMV.df, Conc, Sample)
CONC.DB_by_class_GPMV <- dplyr::bind_rows(DB_agg) #565

save(CONC.DB_by_class_GPMV, file = paste0("./Output_Data/GPMV_CONC.DB_by_Class.Rdata"))
