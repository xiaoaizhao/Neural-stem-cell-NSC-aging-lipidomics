# Format TG in all LC-MS in vitro data, so that it can be used to overlap with Lipidyzer data on correlation analysis

# Overall approach: since TGs reported in Lipidyzer data are annotated with overall carbon and double bond number, we use the same approach ot re-format TGs from LC-MS/MS data
# 0. Reformat LipidIon into ID_string so it can be used to overlap with lipids detected by lipidyzer
# 1. Combine carbon and double bond numbers from all 3 chains
# 2. Check if there are duplicate TG after combining
# 3. For TG with the same total number of carbon and double bond - take the 1 with the highest intensity/concentration

# Note: since annotation output from Lipidyzer is different from LC-MS, we converted lipid annotation from both platforms into a list before overlapping, e.g. c("Cer" ,  "d18:1", "16:0" ) 

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)

source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_Data/Ef_Size_Lipid_Age_qNSC_InVitro_LC.Rdata") # Primary culture #1
load("./Output_Data/Ef_Size_CONC.Lipid_Exp2_all_KO.Rdata") # Primary culture #2, control samples
load("./Output_Data/Exp3_Qui_Age_ES.conc.lipids.Rdata")# Primary culture #3

## ====Exp 1====
E1 <- ether.rename(InVitro.lpd.es.g) %>% 
  rowwise() %>% 
  mutate(., ID_string = list(str_split(LipidIon, "\\(|_|\\)")[[1]][1:length(str_split(LipidIon, "\\(|_|\\)")[[1]])-1])) #373
  
E1.TG.format <- TG.LC.to.overlap(E1, LipidIon)

TGn <- E1.TG.format %>% 
  filter(grepl("^TG", LipidIon)) #51 TG total

DupTG <- E1.TG.format %>% 
  group_by(ID_string_TG) %>% 
  summarise(n = n()) %>% 
  filter(n>1) #need to remove 3 dup TG

TG.clean <- E1.TG.format %>% 
  filter(ID_string_TG %in% DupTG$ID_string_TG) %>% 
  mutate(TG.avg = rowMeans(across(c(Age_mean_Old, Age_mean_Young)))) %>% 
  group_by(ID_string_TG) %>% 
  group_modify(~{
    .x %>% 
      filter(TG.avg == max(.x$TG.avg))
  }) %>% 
  select(-c(ID_string, ID_string_TG.db, ID_string_TG.c, LipidIon1, TG.avg)) %>% 
  rename("ID_string" = "ID_string_TG")

E1.else <- E1 %>% 
  filter(!grepl("^TG", LipidIon)) #322

E1.uniq.TG <- E1.TG.format %>% 
  filter(grepl("TG", LipidIon)) %>% 
  filter(!ID_string_TG %in% DupTG$ID_string_TG) %>%  #45
  select(-c(ID_string, ID_string_TG.db, ID_string_TG.c, LipidIon1)) %>% 
  rename("ID_string" = "ID_string_TG")

E1.fmt <- bind_rows(E1.else, E1.uniq.TG, TG.clean) #370 lipids
save(E1.fmt, file = "./Output_Data/Exp1_EfSz.age.ID_string.reFormat.dupTG.rmved.Rdata")

## ====Exp 2====
Exp2.ctrl.conc <- Exp2.CONC.lpd.es.g.allKO %>% 
  filter(KO == "N") 

E2 <- ether.rename(Exp2.ctrl.conc) %>% 
  rowwise() %>% 
  mutate(., ID_string = list(str_split(LipidIon, "\\(|_|\\)")[[1]][1:length(str_split(LipidIon, "\\(|_|\\)")[[1]])-1])) #612
  
E2.TG.format <- TG.LC.to.overlap(E2, LipidIon)

TGn <- E2.TG.format %>% 
  filter(grepl("^TG", LipidIon)) #43 TG total

DupTG <- E2.TG.format %>% 
  group_by(ID_string_TG) %>% 
  summarise(n = n()) %>% 
  filter(n>1) #need to remove 3 dup TG

TG.clean <- E2.TG.format %>% 
  filter(ID_string_TG %in% DupTG$ID_string_TG) %>% 
  mutate(TG.avg = rowMeans(across(c(Age_mean_Old, Age_mean_Young)))) %>% 
  group_by(ID_string_TG) %>% 
  group_modify(~{
    .x %>% 
      filter(TG.avg == max(.x$TG.avg))
  }) %>% 
  select(-c(ID_string, ID_string_TG.db, ID_string_TG.c, LipidIon1, TG.avg)) %>% 
  rename("ID_string" = "ID_string_TG")

E2.else <- E2 %>% 
  filter(!grepl("^TG", LipidIon)) #569

E2.uniq.TG <- E2.TG.format %>% 
  filter(grepl("TG", LipidIon)) %>% 
  filter(!ID_string_TG %in% DupTG$ID_string_TG) %>%  #45
  select(-c(ID_string, ID_string_TG.db, ID_string_TG.c, LipidIon1)) %>% 
  rename("ID_string" = "ID_string_TG")

E2.fmt <- bind_rows(E2.else, E2.uniq.TG, TG.clean) #609 lipids
save(E2.fmt, file = "./Output_Data/Exp2.ctrl_EfSz.age.ID_string.reFormat.dupTG.rmved.Rdata")

## ====Exp 3====
E3 <- E3.Q.conc.AG.ES %>% 
  rowwise() %>% 
  mutate(., ID_string = list(str_split(LipidIon, "\\(|_|\\)")[[1]][1:length(str_split(LipidIon, "\\(|_|\\)")[[1]])-1]))

E3.TG.format <- TG.LC.to.overlap(E3, LipidIon)

DupTG <- E3.TG.format %>% 
  group_by(ID_string_TG) %>% 
  summarise(n = n()) %>% 
  filter(n>1) #14 unique TG, need to remove 15

TG.clean <- E3.TG.format %>% 
  filter(ID_string_TG %in% DupTG$ID_string_TG) %>% 
  mutate(TG.avg = rowMeans(across(c(Age_mean_Old, Age_mean_Young)))) %>% 
  group_by(ID_string_TG) %>% 
  group_modify(~{
    .x %>% 
      filter(TG.avg == max(.x$TG.avg))
  }) %>% 
  select(-c(ID_string, ID_string_TG.db, ID_string_TG.c, LipidIon1, TG.avg)) %>% 
  rename("ID_string" = "ID_string_TG")

E3.else <- E3 %>% 
  filter(!grepl("^TG", LipidIon)) #268

E3.uniq.TG <- E3.TG.format %>% 
  filter(grepl("TG", LipidIon)) %>% 
  filter(!ID_string_TG %in% DupTG$ID_string_TG) %>%  #57
  select(-c(ID_string, ID_string_TG.db, ID_string_TG.c, LipidIon1)) %>% 
  rename("ID_string" = "ID_string_TG")

E3.fmt <- bind_rows(E3.else, E3.uniq.TG, TG.clean) #339 lipids
save(E3.fmt, file = "./Output_Data/Exp3_EfSz.age.ID_string.reFormat.dupTG.rmved.Rdata")
