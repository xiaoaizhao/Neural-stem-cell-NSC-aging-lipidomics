##This is the first step in double bond composition analysis
##Steps:
##1. Subset by presence or absence of each given double bond, regardless of how many times (side chains) it was identified
##2. Tally total intensity with every given double bond number in each class.
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Pre-processing_functions.R")

####==== All lipids ====
load(file = "./Output_Data/Invivo_Norm_Impt_backtoraw_all130_lipids.Rdata") #this data is in linear concentration

sorted.df <- Invivo.raw_int.invivo %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) #282


##Function to get aggregated intensity for each Class:DB combination
DB_agg <- db.tally(sorted.df, Conc_Int, Sample)
DB_by_class_invivo <- dplyr::bind_rows(DB_agg) #504

save(DB_by_class_invivo, file = paste0("./Output_Data/Invivo_DB_by_Class_qNSC.Rdata"))

####==== Lipids with quantitative standards  ====
#### Side chain analysis on lipids with quantitative standards 
load("./Output_Data/Invivo_Norm_Impt_log2_conc_121_lipids.Rdata")

conc.sorted.df <- 2^Invivo.Impt_norm_conc_all %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) #282

##Function to get aggregated intensity for each Class:DB combination
DB_agg <- db.tally(conc.sorted.df, Conc, Sample)
DB_CONC.by_class_invivo <- dplyr::bind_rows(DB_agg) #504

save(DB_CONC.by_class_invivo, file = paste0("./Output_Data/Invivo_CONC.DB_by_Class_qNSC.Rdata"))
