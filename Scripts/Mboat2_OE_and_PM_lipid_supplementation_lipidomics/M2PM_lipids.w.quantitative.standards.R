##Objective:

# save lipids with quantitative standards from Mboat2 overexpression samples
setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)

library(stringi)

source("./Scripts/Function_scripts/Effect_size_functions.R")

####====Organize experiment M2PM, save quiescent cells only ====####
load("./Output_Data/M2PM_Norm_Impt_log2_conc_275_lipids.Rdata")

M2PM.conc <- 2^M2PM.Impt_norm_conc_all
M2PM.lpd <- M2PM.conc %>% 
  rownames_to_column(var = "LipidIon")

Smp.Key <- read.csv("./Input_Data/March_Sample_list_071123_forR.csv", stringsAsFactors = F)
Smp.Key.e <- Smp.Key %>% 
  select(Sample.Name, Sample_ID)

M2PM.lpd.ctrl <- M2PM.lpd %>% 
  rename_at(vars(matches(Smp.Key.e$Sample_ID)), ~Smp.Key.e$Sample.Name) %>% 
  select(matches("LipidIon|_meth|_EGFP")) %>% 
  mutate(LipidIon = ifelse(grepl("Ch\\+H\\-H2O", LipidIon), "Cholesterol", LipidIon))
M2PM.lpd.ctrl <- ether.rename(M2PM.lpd.ctrl)
save(M2PM.lpd.ctrl, file = "./Output_Data/M2PM.contrl.samples.conc.lipids.Rdata")