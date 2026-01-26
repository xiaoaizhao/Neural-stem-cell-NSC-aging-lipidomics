### Calculate Z score on lipids from 2 in vitro lipidomics
rm(list=ls())
library(RColorBrewer)
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Pre-processing_functions.R")


####Load log2 transformed data from In vitro Experiment #2 dataset####################################################################################
load(file = "./Output_Data/Exp2_MsD.Norm_Impt_log2_conc_278_lipids.Rdata")
df.ctrl <- Exp2.log2.impt.norm.conc.MsD %>%
  rownames_to_column(var = "LipidIon")

##calculate z score for In vitro Experiment #2, control samples only
df.l <- df.ctrl %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Log2Conc") %>% 
  filter(grepl("_N", Sample))

z.ko <- df.l %>%
  group_by(LipidIon) %>% 
  mutate(z = as.vector(scale(Log2Conc))) %>% 
  pivot_wider(-Log2Conc, names_from = Sample, values_from = z) 

save(z.ko, file = "./Output_Data/zscore.lpd.Exp2CTRL.for.cmb.PCA.MsD.frmt.Rdata")
####Load log2 transformed data from In vitro dataset####################################################################################

load("./Output_data/Invitro_MsD.Norm_Impt_log2_conc_404_lipids.Rdata")

Invitro.n <- Invitro.log2.impt.norm.conc.MsD %>% 
  rownames_to_column(var = "LipidIon")

Invitro.l <- Invitro.n %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Log2Conc") %>% 
  filter(grepl("_qNSC-Q", Sample))

z.Invitro <- Invitro.l %>%
  group_by(LipidIon) %>% 
  mutate(z = as.vector(scale(Log2Conc))) %>% 
  pivot_wider(-Log2Conc, names_from = Sample, values_from = z) 


save(z.Invitro, file = "./Output_Data/zscore.lpd.Invitro.for.cmb.PCAMsD.frmt.Rdata")
