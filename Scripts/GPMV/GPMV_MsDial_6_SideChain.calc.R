### Double bond composition analysis

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
source("./Scripts/Function_scripts/Pre-processing_functions.R")

load("./Output_Data/GPMV_MsD.Raw.conc.Norm_Impt_291_lipids.Rdata") #from script #3

## ==== Aggregate lipids with the same number of double bonds for each class, in each sample ====
GPMV.n <- GPMV.impt.norm.raw.conc.lpd.MsD %>% 
  rownames_to_column(var = "LipidIon")

GPMV.rename <- MsD.lpd.format(GPMV.n)

GPMV.SC <- GPMV.rename %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., Class = ifelse(LipidIon == "Cholesterol", "Cholesterol", Class)) %>% 
  mutate(Class = ifelse(Class == "HexCer", "Cer", Class)) %>% 
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) %>% 
  mutate(Sep.SideChain = ifelse(grepl("^Cer|^SM", LipidIon),
                                str_split(SideChain, "/"),
                                str_split(SideChain, "_")))

GPMV.SC.sum <- GPMV.SC %>% 
  unnest(Sep.SideChain) %>%
  mutate(Sep.SC = str_trim(Sep.SideChain)) %>%
  group_by(Class, Sample, Sep.SC) %>% 
  summarise(SumSC = sum(Conc)) %>% 
  filter(!Class == "Cholesterol")

save(GPMV.SC.sum, file = "./Output_Data/GPMV.SC.analysis.Rdata")
