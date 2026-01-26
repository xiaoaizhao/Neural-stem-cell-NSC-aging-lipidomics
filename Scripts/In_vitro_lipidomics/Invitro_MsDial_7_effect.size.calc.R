### In vitro effect size calculation on lipids and double bond composition

setwd(rstudioapi::getActiveProject())
rm(list=ls())
source("./Scripts/Function_scripts/Effect_size_functions.R")
library(tidyverse)
#########################################################################################
load("Output_Data/Invitro_conc.lipid.qNSC.Rdata")
## ====calculate effect size on qNSC from in vitro, all lipids ====####

Invitro.Q <- Invitro.q.conc %>%  
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc")

Invitro.Q.AG.MsD<- Invitro.Q %>%
  mutate(., Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

Invitro.Q.AG.lpd.es <- es.g.func(Invitro.Q.AG.MsD, LipidIon, Age, Conc, Sample)

## ==== Append RT to Effect size matrix ==== 

Invitro<- read.csv("./Input_Data/Exp3_Validated_Reintegration_121125_color.code.csv", stringsAsFactors = F)
load("./Output_Data/Invitro.MS-Dial_IS.all.Rdata")

Invitro.lpd.check <- Invitro %>% 
  rowwise() %>% 
  mutate(Class = str_split(Metabolite.name, " ")[[1]][1]) %>% 
  mutate(Class = ifelse(Class == "ST", "Cholesterol", Class)) %>% 
  group_by(Metabolite.name) %>% 
  group_modify(~{
    .x %>% 
      mutate(Ion.check = ifelse(Adduct.type == IS.invitro.export.flt$Adduct.type[IS.invitro.export.flt$Class == Class], "T", "F"))  
  }) %>% 
  relocate(c(Class, Ion.check), .after = "Metabolite.name") %>% 
  filter(Ion.check == "T")  %>% ##all 441 lipids check out
  select(Metabolite.name, Average.Rt.min.) %>% 
  rename("LipidIon" = "Metabolite.name")

Invitro.Q.lpd.ES.w.RT <- left_join(Invitro.Q.AG.lpd.es, Invitro.lpd.check, by = "LipidIon") %>% 
  relocate(Average.Rt.min., .after = "LipidIon")

save(Invitro.Q.lpd.ES.w.RT, file = paste0("./Output_data/Lpd.MsD.Invitro_Qui_Age_ES.w.RT.Rdata"))
## ====calculate effect size on qNSC from Experiment #3, Double bond composition with quantitative standards====####
load("./Output_Data/Invitro.Qui_CONC.DB_PCT_by_Class.Rdata")
Qui_CONC.DB.invitro <- MsD.Qui_CONC.DB.invitro %>%
  mutate(., Cla_DB = paste0(Class, DB_num))

Invitro.Qui.MsD.CONC.DB.es.g <- es.g.func(Qui_CONC.DB.invitro, Cla_DB, Age, DB_Pct, Sample)
save(Invitro.Qui.MsD.CONC.DB.es.g, file = "./Output_Data/DBPct.MsD.Invitro_Qui_Age_ES.Rdata")

