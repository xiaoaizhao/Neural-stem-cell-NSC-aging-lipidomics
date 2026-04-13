### Side chain composition analysis

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
source("./Scripts/Function_scripts/Pre-processing_functions.R")

load("./Output_Data/Invivo.MsD.29.lipid.raw.conc.Rdata") 

## ==== Aggregate lipids with the same number of double bonds for each class, in each sample ====
Invivo.rename <- MsD.lpd.format(Invivo.raw.conc.29.good.peak)
Invivo.SC <- Invivo.rename %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., Class = ifelse(LipidIon == "Cholesterol", "Cholesterol", Class)) %>% 
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) %>% 
  mutate(Sep.SideChain = ifelse(grepl("^Cer|^SM", LipidIon),
                                str_split(SideChain, "/"),
                                str_split(SideChain, "_")))

Invivo.SC.sum <- Invivo.SC %>% 
  unnest(Sep.SideChain) %>%
  mutate(Sep.SC = str_trim(Sep.SideChain)) %>%
  group_by(Class, Sample, Sep.SC) %>% 
  summarise(SumSC = sum(Conc)) %>% 
  filter(!Class == "Cholesterol")

save(Invivo.SC.sum, file = "./Output_Data/Invivo.SC.analysis.Rdata")

