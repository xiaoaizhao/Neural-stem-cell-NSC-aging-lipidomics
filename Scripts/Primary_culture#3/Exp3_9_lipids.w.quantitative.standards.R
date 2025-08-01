##Objective:

# save lipids with quantitative standards from Experiment 3 data, separate quiescent cells only, as well as quiescent and activated cells together.

rm(list=ls())
library(tidyverse)

library(stringi)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")

####====Organize experiment 3, save quiescent cells only ====####
load("./Output_Data/Exp3_Norm_Impt_log2_conc_356_lipids.Rdata")

B3.LpdSearch <- 2^Exp3.Impt_norm_conc_all

B3.LS.Cla<- B3.LpdSearch%>% 
  rownames_to_column(var = "LipidIon") 

Smp.Key <- read.csv("./Input_Data/Batch3_sample_key_forR.csv", stringsAsFactors = F)

Smp.Key.e <- Smp.Key %>% 
  filter(!grepl("PBS", Sample.Name))

B3.Q.Cla<- B3.LS.Cla %>% 
  rename_at(vars(matches(Smp.Key.e$Sample_ID)), ~Smp.Key.e$Sample.Name) %>% 
  select(-"O8_aNSC-A") %>% 
  select(matches("_qNSC-Q|LipidIon"))

Exp3.conc.Q <- B3.Q.Cla %>% 
  filter(!grepl("Cholesterol\\+H\\-H2O_23_Positive|Ch\\+H\\-H2O", LipidIon)) 

save(Exp3.conc.Q, file = "./Output_Data/Exp3_qNSC_quant.lipids.Rdata")


####====Save activated and quiescent cells ====####
Exp3.conc.QA <- B3.LS.Cla %>% 
  rename_at(vars(matches(Smp.Key.e$Sample_ID)), ~Smp.Key.e$Sample.Name) %>% 
  select(-"O8_aNSC-A") %>% 
  filter(!grepl("Cholesterol\\+H\\-H2O_23_Positive|Ch\\+H\\-H2O", LipidIon)) 

save(Exp3.conc.QA, file = "./Output_Data/Exp3_qNSC.aNSC_quant.lipids.Rdata")


####====All lipids in all samples ====####
load(file = "./Output_data/Exp3_backtoRAW_Norm_Impt_log2_conc_356_noconc_60_lipids.Rdata")

Smp.Key <- read.csv("~/Dropbox/Stanford/Lipidomics/2023_Feb_lipidomics/Batch3_sample_key_forR.csv", stringsAsFactors = F)
Smp.Key.e <- Smp.Key %>% 
  select(Sample.Name, Sample_ID) %>% 
  filter(!grepl("PBS", Sample_ID))

Exp3.all.lpd.all.cell <- raw_int.Exp3 %>% 
  rename_at(vars(matches(Smp.Key.e$Sample_ID)), ~Smp.Key.e$Sample.Name) %>% 
  select(-`O8_aNSC-A`) %>% 
  rownames_to_column(var = "LipidIon") %>% 
  filter(!grepl("Cholesterol\\+H\\-H2O_23_Positive|Ch\\+H\\-H2O", LipidIon))

save(Exp3.all.lpd.all.cell, file = "./Output_Data/Exp3_all.lipid.all.cells.Rdata")