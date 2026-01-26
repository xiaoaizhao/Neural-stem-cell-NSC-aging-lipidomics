### In vivo effect size calculation on lipids and double bond composition

setwd(rstudioapi::getActiveProject())
rm(list=ls())
source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")
library(tidyverse)
#########################################################################################

## ====calculate effect size on qNSC from Experiment #3, all lipids ====####
load("Output_Data/Invivo.MsD.29.lipid.raw.conc.Rdata")
Invivo <- Invivo.raw.conc.29.good.peak %>%  
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc")


Invivo.AG.MsD<- Invivo %>%
  mutate(., Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

Invivo.AG.MsD.ES.lpd<- es.g.func(Invivo.AG.MsD, LipidIon, Age, Conc, Sample)

## ==== Append RT to Effect size matrix ==== 
load("./Output_Data/Invivo_IS_w_conc.Rdata")
Invivo.df <- read.csv("./Input_Data/Invivo_sorted_NSCs_integration_graded_120325.csv", stringsAsFactors = F)

df.fix <- Invivo.df %>% 
  mutate(Metabolite.name = ifelse(Alignment.ID == 3476, "ST 27:1;O", Metabolite.name)) %>% 
  mutate(Metabolite = ifelse(Alignment.ID == 3476, "ST 27:1;O", Metabolite)) 

Invivo.lpd.check <- df.fix %>% 
  rowwise() %>% 
  mutate(Class = str_split(Metabolite.name, " ")[[1]][1]) %>%
  mutate(Class = ifelse(Class == "ST", "Cholesterol", Class)) %>% 
  relocate(Class, .after = Metabolite.name) %>% 
  group_by(Metabolite.name) %>% 
  group_modify(~{
    .x %>% 
      mutate(Ion.check = ifelse(Adduct.type == Invivo.IS.Conc$Adduct.type[Invivo.IS.Conc$Class == Class], "T", "F"))  
  }) %>% 
  relocate(c(Class, Ion.check), .after = "Metabolite.name") %>% 
  filter(Ion.check == "T") %>% 
  ungroup() %>% 
  select(Metabolite, Average.Rt.min.) %>% 
  rename("LipidIon" = "Metabolite")

Invivo.lpd.ES.w.RT <- left_join(Invivo.AG.MsD.ES.lpd, Invivo.lpd.check, by = "LipidIon") %>% 
  relocate(Average.Rt.min., .after = "LipidIon")

save(Invivo.lpd.ES.w.RT, file = paste0("./Output_data/Lpd.MsD.Invivo_Age_ES.w.RT.Rdata"))

## ====calculate effect size on qNSC from Experiment #3, Double bond composition with quantitative standards====####
load("./Output_Data/MsD.Invivo_CONC.DB_PCT_by_Class.Rdata")
DB.pct.invivo <- MsD.DB.pct.invivo %>%
  mutate(., Cla_DB = paste0(Class, DB_num))

Invivo.DB.AG.ES<- es.g.func(DB.pct.invivo, Cla_DB, Age, DB_Pct, Sample)

save(Invivo.DB.AG.ES, file = "./Output_Data/DBPct.MsD.Invivo_Age_ES.Rdata")

