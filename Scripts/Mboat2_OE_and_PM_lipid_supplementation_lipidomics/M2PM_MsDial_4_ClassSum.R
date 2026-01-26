### Class sum Mboat2 overexpression and plasma membrane lipid supplementation
setwd(rstudioapi::getActiveProject())
library(tidyverse)
rm(list=ls())

load("./Output_Data/M2PM_backtoRAW_MsD_Norm_Impt_397_lipids.Rdata")

M2PM.class.sum <- M2PM.MsD.lpd.conc.all %>% 
  rownames_to_column(var = "LipidIon") %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Concentration") %>% 
  rowwise() %>% 
  mutate(Class = str_split(LipidIon, " ")[[1]][1]) %>% 
  mutate(Class = ifelse(Class == "ST", "Cholesterol", Class)) %>% 
  mutate(Class = ifelse(Class == "HexCer", "Cer", Class)) %>% 
  group_by(Class, Sample) %>% 
  summarise(ClassSum = sum(Concentration))

save(M2PM.class.sum, file = "./Output_Data/M2PM.MsD.ClassSum.Rdata")

M2PM.total.conc <- M2PM.class.sum %>% 
  group_by(Sample) %>% 
  summarise(TotalConc = sum(ClassSum))

save(M2PM.total.conc, file = "./Output_Data/M2PM.MsD.TotalLipid.Rdata")
