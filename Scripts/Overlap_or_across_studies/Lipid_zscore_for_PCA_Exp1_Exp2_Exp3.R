## PCA on all 3 datasets together, Exp1, Exp2 control and Batch 3
## Re-format ether lipid in Exp1 and Exp2 data
## Use the endogenous cholesterol peak eluded at ~10min from Xcalibur for all downstream analysis - remove cholesterol annotated by LipidSearch and a second cholesterol annotated by Xcalibur that has an RT of ~23mins
## Identify commonly detected lipids across all datasets (95 common lipids)
## PCA analysis

rm(list=ls())
library(RColorBrewer)
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")
####Load log2 transformed data from 1st Primary NSC culture data####################################################################################
load(file = "./Output_data/Spike-in_norm_MedNorm_all_373_lipid.Rdata") #this data is already log2 transformed

##calculate z score for Primary Culture #1 dataset
Exp1 <- Impt_norm_373all %>% 
  rownames_to_column(var = "LipidIon") 

Exp1.n <- ether.rename(Exp1) 
  

z.lc <- Exp1.n %>%
  pivot_longer(-LipidIon, names_to = "Samples", values_to = "Int") %>%
  group_by(LipidIon) %>% 
  mutate(z = as.vector(scale(Int))) %>% 
  pivot_wider(-Int, names_from = Samples, values_from = z) 

save(z.lc, file = "./Output_Data/zscore.lpd.Exp1.for.cmb.PCA.Rdata")

####Load log2 transformed data from Primary Culture #2 dataset####################################################################################
load(file = "./Output_Data/Exp2_Norm_Impt_log2_all693_lipids.Rdata")
df.ctrl <- Exp2_Impt_norm_conc_no_conc_all %>%
  select(., contains("_N"))

##calculate z score for Primary Culture #2 dataset, control samples only
Exp2 <- df.ctrl %>% 
  rownames_to_column(var = "LipidIon")
Exp2.n <- ether.rename(Exp2) 
  
z.ko <- Exp2.n %>%
  pivot_longer(-LipidIon, names_to = "Samples", values_to = "Int") %>%
  group_by(LipidIon) %>% 
  mutate(z = as.vector(scale(Int))) %>% 
  pivot_wider(-Int, names_from = Samples, values_from = z) 

save(z.ko, file = "./Output_Data/zscore.lpd.Exp2CTRL.for.cmb.PCA.Rdata")
####Load log2 transformed data from Primary culture #3 dataset####################################################################################
load("./Output_data/Exp3_Norm_Impt_log2_conc_356_noconc_60_lipids.Rdata")

Smp.Key <- read.csv("./Input_Data/Batch3_sample_key_forR.csv", stringsAsFactors = F)
Smp.Key.e <- Smp.Key %>% 
  select(Sample.Name, Sample_ID) %>% 
  filter(!grepl("PBS", Sample_ID))

## Annotate sample ID with sample name
Exp3.df <- Exp3.Impt_norm_conc_no_conc_all %>% 
  rename_at(vars(matches(Smp.Key.e$Sample_ID)), ~Smp.Key.e$Sample.Name)

Exp3.no.O8A <- Exp3.df %>% 
  select(-`O8_aNSC-A`) 

##calculate z score for Batch#3 dataset, control samples only
Exp3 <- Exp3.no.O8A %>% 
  rownames_to_column(var = "LipidIon") %>% 
  filter(!grepl("Cholesterol\\+H\\-H2O_23_Positive|Ch\\+H\\-H2O", LipidIon))

Exp3.n <- ether.rename(Exp3)

z.Exp3 <- Exp3.n %>%
  pivot_longer(-LipidIon, names_to = "Samples", values_to = "Int") %>%
  group_by(LipidIon) %>% 
  mutate(z = as.vector(scale(Int))) %>% 
  pivot_wider(-Int, names_from = Samples, values_from = z) 


save(z.Exp3, file = "./Output_Data/zscore.lpd.Exp3.for.cmb.PCA.Rdata")
