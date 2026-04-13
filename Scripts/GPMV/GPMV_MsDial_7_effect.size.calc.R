## In vivo sorted cells effect size calculation

setwd(rstudioapi::getActiveProject())
rm(list=ls())
source("./Scripts/Function_scripts/Effect_size_functions.R")
library(tidyverse)
#########################################################################################

## ====calculate effect size on GPMV, all lipids ====####
load("Output_Data/GPMV_MsD.Raw.conc.Norm_Impt_291_lipids.Rdata")
GPMV <- GPMV.impt.norm.raw.conc.lpd.MsD %>%  
  rownames_to_column(var = "LipidIon") %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc")


GPMV.AG.MsD<- GPMV %>%
  mutate(., Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

GPMV.AG.MsD.ES.lpd <- es.g.func(GPMV.AG.MsD, LipidIon, Age, Conc, Sample)

## ==== Append RT to Effect size matrix ==== 
load("./Output_Data/GPMV.IS.w.conc.Rdata")
GPMV.df <- read.csv("./Input_Data/GPMV_validated_120525_Full_integration_peak_graded.csv")

GPMV.lpd.check <- GPMV.df %>% 
  filter(!Poor_peak_shape == "Yes") %>% 
  rowwise() %>%
  mutate(Class = str_split(Metabolite.name, " ")[[1]][1]) %>%
  mutate(Class = ifelse(Class == "ST", "Cholesterol", Class)) %>%
  mutate(Class = ifelse(Class == "HexCer", "Cer", Class)) %>%
  group_by(Metabolite.name) %>%
  group_modify(~{
    .x %>%
      mutate(Ion.check = ifelse(Adduct.type == GPMV.IS.w.conc$Adduct.type[GPMV.IS.w.conc$Class == Class], "T", "F"))
  }) %>%
  relocate(c(Class, Ion.check), .after = "Metabolite.name") %>% 
  filter(Ion.check == "T") %>% 
  ungroup() %>% 
  select(Metabolite.name, Average.Rt.min.) %>% 
  rename("LipidIon" = "Metabolite.name")

GPMV.lpd.ES.w.RT.291.good.peak <- left_join(GPMV.AG.MsD.ES.lpd, GPMV.lpd.check, by = "LipidIon") %>% 
  relocate(Average.Rt.min., .after = "LipidIon")


save(GPMV.lpd.ES.w.RT.291.good.peak, file = paste0("./Output_data/Lpd.MsD.GPMV_Age_ES.291lpd.Rdata"))

## ====calculate effect size on qNSC from GPMV, side chain analysis====####
load("./Output_Data/GPMV.SC.analysis.Rdata")

SC.GPMV <- GPMV.SC.sum %>% 
  mutate(Cla_SC = paste0(Class, "(", Sep.SC, ")")) %>% 
  mutate(Age = ifelse(grepl("^O", Sample), "Old", "Young")) 

GPMV.SC.es.g <- es.g.func(SC.GPMV, Cla_SC, Age, SumSC, Sample)
save(GPMV.SC.es.g, file = "./Output_Data/SC.abundance.GPMV_Age_ES.Rdata")

