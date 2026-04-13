## In vivo sorted cells effect size calculation

setwd(rstudioapi::getActiveProject())
rm(list=ls())
source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")
library(tidyverse)
#########################################################################################

## ====calculate effect size on qNSC in vivo, all lipids ====####
load("./Output_Data/Invivo.MsD.29.lipid.raw.conc.Rdata")
Invivo <- Invivo.raw.conc.29.good.peak %>%  
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc")


Invivo.AG.MsD<- Invivo %>%
  mutate(., Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

Invivo.AG.MsD.ES.lpd <- es.g.func(Invivo.AG.MsD, LipidIon, Age, Conc, Sample)

## ==== Append RT to Effect size matrix ==== 
load("./Output_Data/Invivo_IS_w_conc.Rdata")
data_book1<-loadWorkbook("./Input_Data/Invivo_sorted_NSCs_integration_cleanup_112325.xlsx")
data<-list()

for (i in sheets(data_book1)){
  data[[i]]<-data.frame(read.xlsx("./Input_Data/Invivo_sorted_NSCs_integration_cleanup_112325.xlsx",sheet=i))
}

Invivo.df <- bind_rows(data)

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
  filter(Ion.check == "T") %>% ##all 51 lipids check out
  ungroup() %>% 
  select(Metabolite, Average.Rt.min.) %>% 
  rename("LipidIon" = "Metabolite")

Invivo.lpd.ES.w.RT <- left_join(Invivo.AG.MsD.ES.lpd, Invivo.lpd.check, by = "LipidIon") %>% 
  relocate(Average.Rt.min., .after = "LipidIon")

save(Invivo.lpd.ES.w.RT, file = paste0("./Output_Data/Lpd.MsD.Invivo_Age_ES.w.RT.Rdata"))


## ====calculate effect size on qNSC from In vivo, side chain analysis====####
load("./Output_Data/Invivo.SC.analysis.Rdata")

SC.invivo <- Invivo.SC.sum %>% 
  mutate(Cla_SC = paste0(Class, "(", Sep.SC, ")")) %>% 
  mutate(Age = ifelse(grepl("^O", Sample), "Old", "Young")) 

Invivo.SC.es.g <- es.g.func(SC.invivo, Cla_SC, Age, SumSC, Sample)
save(Invivo.SC.es.g, file = "./Output_Data/SC.abundance.Invivo_Age_ES.Rdata")
