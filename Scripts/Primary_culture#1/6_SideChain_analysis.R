##This is the first step in double bond composition analysis
##Steps:
##1. Subset by presence or absence of each given double bond, regardless of how many times (side chains) it was identified
##2. Tally total intensity with every given double bond number in each class.

rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load(file = "./Output_data/Spike-in_norm_Mednorm_all_373_lipids_back_to_raw_int.Rdata")
WC <- raw_int %>%
  mutate(., LipidIon = rownames(raw_int))
WC<- WC %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Int") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) #8952/373=24, include both quiescent and activated samples

##Function to get aggregated intensity for each Class:DB combination####
DB_agg <- db.tally(WC, Int, Sample)
DB_by_class_2017 <- dplyr::bind_rows(DB_agg) #1608

save(DB_by_class_2017, file = paste0("./Output_Data/Act_Qui_DB_By_Class_2017LC.Rdata"))

##Quiescent samples only
Qui_DB_By_Cla_2017 <- DB_by_class_2017 %>% #804
  filter(., grepl("Qui", Sample, ignore.case = F))

save(Qui_DB_By_Cla_2017, file = paste0("./Output_data/Qui_DB_by_Class_2017LC.Rdata"))

