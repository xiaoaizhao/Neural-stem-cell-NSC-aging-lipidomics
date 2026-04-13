### Side Chain composition analysis

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
source("./Scripts/Function_scripts/Pre-processing_functions.R")

load("./Output_Data/M2PM_backtoRAW_MsD_Norm_Impt_397_lipids.Rdata") #from script #6

## ==== Aggregate lipids with the same number of double bonds for each class, in each sample ====
M2PM.n <- M2PM.MsD.lpd.conc.all %>% 
  rownames_to_column(var = "LipidIon")

M2PM.rename <- MsD.lpd.format(M2PM.n)

M2PM.SC <- M2PM.rename %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., Class = ifelse(LipidIon == "Cholesterol", "Cholesterol", Class)) %>% 
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) %>% 
  mutate(Sep.SideChain = ifelse(grepl("^Cer|^SM|^HexCer", LipidIon),
                                str_split(SideChain, "/"),
                                str_split(SideChain, "_")))

M2PM.SC.sum <- M2PM.SC %>% 
  unnest(Sep.SideChain) %>%
  mutate(Sep.SC = str_trim(Sep.SideChain)) %>%
  group_by(Class, Sample, Sep.SC) %>% 
  summarise(SumSC = sum(Conc))

save(M2PM.SC.sum, file = "./Output_Data/M2PM.SC.analysis.Rdata")


### Separate Mboat2 OE and PM lipid supplementation samples
M2OE.SC.sum <- M2PM.SC.sum %>% 
  filter(grepl("_EGFP|_Mb2_OE", Sample)) #3184
save(M2OE.SC.sum, file = "./Output_Data/Mboat2OE.SC.analysis.Rdata")
PM.sup.SC.sum <- M2PM.SC.sum %>% 
  filter(grepl("_meth|lpd", Sample)) #4776
save(PM.sup.SC.sum, file = "./Output_Data/PM.lipid.supp.SC.analysis.Rdata")