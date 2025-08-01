## Prerequisite: Need to run `M2PM_lipids.w.quantitative.standards.R` first
## Quantification on the class concentration of control samples from Mboat2 overexpression samples
rm(list = ls())
library(tidyverse)
library(ggbeeswarm)
library(ggpubr)
source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_Data/M2PM.contrl.samples.conc.lipids.Rdata")

M2PM.cla <- M2PM.lpd.ctrl %>% 
  rowwise() %>% 
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>% 
  mutate(., Class = ifelse(grepl("Cholesterol", LipidIon), "Cholesterol", Class))


M2PM.ctrl.cla.sum <- M2PM.cla %>% 
  pivot_longer(-c(LipidIon, Class), names_to = "Samples", values_to = "Concentration") %>% 
  group_by(Samples, Class) %>% 
  summarise(., ClassSum = sum(Concentration)) %>% 
  mutate(., Age = ifelse(grepl("Y", Samples), "Young", "Old"))


save(M2PM.ctrl.cla.sum, file = "./Output_Data/M2PM.CTRL.ClassConcSum.Rdata")
