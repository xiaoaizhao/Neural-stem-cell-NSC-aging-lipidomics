
# Concentration class sum difference between young and old quiescent cells from GPMV dataset
# Only included lipis with quantitative standards
# Remove mis-annotated LPC from the dataset
rm(list = ls())
library(tidyverse)
library(ggbeeswarm)
library(ggpubr)

source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_Data/Invivo_Norm_Impt_RAW.conc_121_lipids.Rdata")

## list of misannotated LPC from manual validation by tandem mass spectrometry
lpc.to.rmv <- c("LPC(14:0)", "LPC(16:0e)", "LPC(16:1e)", "LPC(16:2e)", "LPC(17:1)", "LPC(18:1e)", "LPC(18:3e)", "LPC(20:1)", "LPC(20:5)")

Invivo.Cla.org <- Invivo.Impt.norm.conc.raw %>%
  rowwise() %>% 
  mutate(LipidID = substr(LipidIon, 1, str_locate(LipidIon, "\\)"))) %>% 
  filter(!LipidID %in% lpc.to.rmv) %>% 
  select(-LipidID)

Invivo.Cla.org.l <- Invivo.Cla.org %>% 
  pivot_longer(-c(LipidIon, Class), names_to = "Samples", values_to = "Concentration") %>% 
  group_by(Samples, Class) %>% 
  summarise(., ClassSum = sum(Concentration)) %>% 
  mutate(., Age = ifelse(grepl("Y", Samples), "Young", "Old"))

save(Invivo.Cla.org.l, file = paste0("./Output_Data/Invivo.ClassConcSum.Rdata"))

