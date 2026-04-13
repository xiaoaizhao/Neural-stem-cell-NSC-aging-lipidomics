### In vitro #2 effect size calculation on lipids and double bond composition

setwd(rstudioapi::getActiveProject())
rm(list=ls())
source("./Scripts/Function_scripts/Effect_size_functions.R")
library(tidyverse)
#########################################################################################

## ====calculate effect size on qNSC from In vitro #2, all lipids ====####
load("Output_Data/Exp2_278_conc_lipid.MSD.format.Rdata")
E2 <- E2MsD.frmt %>%  
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc")


E2.KO.MsD<- E2 %>%
  mutate(., Age = ifelse(grepl("^Y", Sample), "Young", "Old")) %>% 
  mutate(., KO = substr(Sample, nchar(Sample), nchar(Sample)))

##effect size on age difference at individual lipid level within each KO condition
KO.list <- unique(E2.KO.MsD$KO)

EfSize.df = list()
for (KO.name in KO.list) {
  KO.df<- E2.KO.MsD %>%
    filter(., KO == KO.name)
  EfSize.df[[KO.name]] <- es.g.func(KO.df, LipidIon, Age, Conc, Sample) %>%
    mutate(., KO = KO.name)
}

Exp2.lpd.es.g.MsD.allKO <- bind_rows(EfSize.df)

## ==== Append RT to Effect size matrix ==== 
Exp2 <- read.csv("./Input_Data/Exp2_Reintegration_values_121725.peak.color.coded.csv", stringsAsFactors = F)
load("./Output_Data/Exp2.MS-Dial_IS.all.Rdata")

E2.lpd.check <- Exp2 %>% 
  rowwise() %>% 
  mutate(Class = str_split(Metabolite.name, " ")[[1]][1]) %>% 
  group_by(Metabolite.name) %>% 
  group_modify(~{
    .x %>% 
      mutate(Ion.check = ifelse(Adduct.type == IS2.export.flt$Adduct.type[IS2.export.flt$Class == Class], "T", "F"))  
  }) %>% 
  relocate(c(Class, Ion.check), .after = "Metabolite.name") %>% 
  filter(Ion.check == "T") %>% ##all 287 lipids checked out
  select(Metabolite.name, Average.Rt.min.) %>% 
  rename("LipidIon" = "Metabolite.name")

Exp2.lpd.es.g.MsD.allKO.wRT <- left_join(Exp2.lpd.es.g.MsD.allKO, E2.lpd.check, by = "LipidIon") %>% 
  relocate(Average.Rt.min., .after = "LipidIon")

save(Exp2.lpd.es.g.MsD.allKO.wRT, file = "./Output_Data/Lpd.MsD.Age.ES_Exp2_all_KO.w.RT.Rdata")


## ====calculate effect size on all KO from In vitro Experiment #2, side chain abundance ====####
load("./Output_Data/E2.SC.analysis.Rdata")
Exp2_SC <- E2.SumSC %>%
  mutate(Cla_SC = paste0(Class, "(", Sep.SC, ")")) %>% 
  mutate(., KO = substr(Samples, nchar(Samples), nchar(Samples))) %>% 
  mutate(Age = ifelse(grepl("^O", Samples), "Old", "Young")) 

##effect size on age difference at double bond composition level within each KO condition
KO.list <- unique(Exp2_SC$KO)

EfSize.df.SC = list()
for (KO.name in KO.list) {
  KO.df<- Exp2_SC %>%
    filter(., KO == KO.name)
  EfSize.df.SC[[KO.name]] <- es.g.func(KO.df, Cla_SC, Age, SumSC, Samples) %>%
    mutate(., KO = KO.name)
}

Exp2.SC.es.g.allKO <- bind_rows(EfSize.df.SC)
save(Exp2.SC.es.g.allKO, file = "./Output_Data/SC.abundance.MsD.Age.ES_Exp2_all_KO.Rdata")
