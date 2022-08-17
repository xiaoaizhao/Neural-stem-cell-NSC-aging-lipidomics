
##This is the first step in double bond composition analysis
##Steps:
##1. Subset by presence or absence of each given double bond, regardless of how many times (side chains) it was identified
##2. Tally total intensity with every given double bond number in each class.
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Pre-processing_functions.R")

###Import Exp2 lipid normalized with imputation, linear concentration##################################################################################
###Import Exp2 lipid normalized with imputation, linear concentration##################################################################################
load(file = "./Output_Data/Exp2_Norm_Impt_backtoraw_all693_lipids.Rdata")

KO.df <- raw_int.exp2 %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) #7904/494 = 16

##Function to get aggregated intensity for each Class:DB combination####
DB_agg <- db.tally(KO.df, Conc_Int, Sample)
Exp2_DB_by_class <- dplyr::bind_rows(DB_agg) #3840

save(Exp2_DB_by_class, file = paste0("./Output_data/Exp2_DB_by_Class.Rdata"))

