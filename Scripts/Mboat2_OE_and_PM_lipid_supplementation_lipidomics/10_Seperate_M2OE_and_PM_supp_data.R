# Separate DB PCT and raw lipid data into samples of _Mboat2_ overexpression and young plasma membrane lipid supplementation samples
## --------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(RColorBrewer)
library(rstatix)

setwd(rstudioapi::getActiveProject())
## -------------------------------------------------------------------------------------------------------------------
load('Output_Data/M2PM_all_sample_DB_PCT_by_Class.Rdata') # effect size matrix on double bond composition in Batch 2 data with Mboat2 OE and GPMV supplementation
load("./Output_Data/M2PM_backtoRAW_Norm_Impt_log2_conc_275_noconc_48_lipids.Rdata") # effect size matrix on individual lipid in Primary Culture #2 Exp with KO


## Separate DB PCT of 2 studies, Mboat2 OE and GPMV supplementation
## -------------------------------------------------------------------------------------------------------------------
Smp.Key <- read.csv("./Input_Data/March_Sample_list_071123_forR.csv", stringsAsFactors = F)
Smp.Key.e <- Smp.Key %>% 
  select(Sample.Name, Sample_ID)  %>% 
  rename("Sample" = "Sample_ID")

M2PM.db <- left_join(M2PM.DB.pct, Smp.Key.e, by = "Sample")

## Name check
n=1085
M2PM.db$Sample[n] == Smp.Key.e$Sample[Smp.Key.e$Sample.Name == M2PM.db$Sample.Name[n]]

M2OE.DB <- M2PM.db %>% 
  filter(grepl("_EGFP|_Mb2_OE", Sample.Name)) #1136/16 = 71

PM.sup.DB <- M2PM.db %>% 
  filter(grepl("_meth|lpd", Sample.Name)) #1704/24 = 71

save(M2OE.DB, file = "./Output_Data/Mboat2_OE_DB_PCT.Rdata")
save(PM.sup.DB, file = "./Output_Data/PM_supp_DB_PCT.Rdata")

## Separate individual lipid data matrix of 2 studies, Mboat2 OE and GPMV supplementation
## -------------------------------------------------------------------------------------------------------------------

M2PM.lipid <- M2PM.raw.all %>% 
  rename_at(vars(matches(Smp.Key.e$Sample)), ~Smp.Key.e$Sample.Name)

## Name check
n=17
colnames(M2PM.raw.all)[n] == Smp.Key.e$Sample[Smp.Key.e$Sample.Name == colnames(M2PM.lipid)[n]]

M2OE.Lipid <- M2PM.lipid %>% 
  select(matches("_EGFP|_Mb2_OE")) # 16 samples

PM.sup.Lipid <- M2PM.lipid %>% 
  select(matches("_meth|lpd")) # 24 samples

save(M2OE.Lipid, file = "./Output_Data/Mboat2_OE_LIPID.Rdata")
save(PM.sup.Lipid, file = "./Output_Data/PM_supp_LIPID.Rdata")
