## Primary culture #3 effect size calculation

setwd(rstudioapi::getActiveProject())
rm(list=ls())
source("./Scripts/Function_scripts/Effect_size_functions.R")
library(tidyverse)
#########################################################################################

## ====calculate effect size on qNSC from Experiment #3, all lipids ====####
load("Output_Data/Exp3_all.lipid.all.cells.Rdata")
E3.Q <- ether.rename(Exp3.all.lpd.all.cell) %>% 
  select(matches("LipidIon|qNSC-Q")) %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int")


E3.Q.AG<- E3.Q %>%
  mutate(., Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

E3.Q.AG.ES <- es.g.func(E3.Q.AG, LipidIon, Age, Conc_Int, Sample)

save(E3.Q.AG.ES, file = paste0("./Output_data/Exp3_Qui_Age_ES.Rdata"))

## ====calculate effect size on qNSC from Experiment #3, lipids with quantitative standards====####
load("./Output_Data/Exp3_qNSC_quant.lipids.Rdata")
E3.Q.conc <- ether.rename(Exp3.conc.Q) %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc")


E3.Q.conc.AG<- E3.Q.conc %>%
  mutate(., Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

E3.Q.conc.AG.ES <- es.g.func(E3.Q.conc.AG, LipidIon, Age, Conc, Sample)

save(E3.Q.conc.AG.ES, file = paste0("./Output_data/Exp3_Qui_Age_ES.conc.lipids.Rdata"))

## ====calculate effect size between qNSCs and aNSCs from Experiment #3, lipids with quantitative standards====####
load("./Output_Data/Exp3_all.lipid.all.cells.Rdata")
B3.all <- ether.rename(Exp3.all.lpd.all.cell) %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int")

Cell.type.Invitro <- B3.all %>%
  mutate(., CellType = ifelse(grepl("aNSC-A", Sample), "Activated", "Quiescent"))

Exp3.CellType.ES <- es.g.func.celltype(Cell.type.Invitro, LipidIon, CellType, Conc_Int, Sample) 

save(Exp3.CellType.ES, file = paste0("./Output_Data/Exp3_Cell_Type_ES.Rdata"))

## ====calculate effect size on qNSC from Experiment #3, Double bond composition with quantitative standards====####
load("./Output_Data/Exp3.Qui_CONC.DB_PCT_by_Class.Rdata")
Qui_CONC.DB.E3 <- Qui_CONC.DB.E3 %>%
  mutate(., Cla_DB = paste0(Class, DB_num))

Exp3.Qui.CONC.DB.es.g <- es.g.func(Qui_CONC.DB.E3, Cla_DB, Age, DB_Pct, Sample)

save(Exp3.Qui.CONC.DB.es.g, file = "./Output_Data/Ef_Size_CONC.DB_pct_Exp3_Qui_aging.Rdata")

## ====calculate effect size on qNSC from Experiment #3, All Double bond composition====####
load("./Output_Data/Exp3.Qui_all.DB_PCT_by_Class.Rdata")
Qui.DB.E3 <- Qui.DB.E3 %>%
  mutate(., Cla_DB = paste0(Class, DB_num))

Exp3.Qui.DB.es.g <- es.g.func(Qui.DB.E3, Cla_DB, Age, DB_Pct, Sample)

save(Exp3.Qui.DB.es.g, file = "./Output_Data/Ef_Size_ALL.DB_pct_Exp3_Qui_aging.Rdata")
