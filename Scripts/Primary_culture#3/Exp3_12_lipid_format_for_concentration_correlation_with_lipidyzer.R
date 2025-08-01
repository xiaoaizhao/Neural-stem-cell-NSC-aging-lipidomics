
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(stringi)

source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_Data/Exp3_qNSC.aNSC_quant.lipids.Rdata")

E3 <- Exp3.conc.QA %>% 
  rowwise() %>% 
  mutate(., ID_string = list(str_split(LipidIon, "\\(|_|\\)")[[1]][1:length(str_split(LipidIon, "\\(|_|\\)")[[1]])-1]))

E3.TG.format <- TG.LC.to.overlap(E3, LipidIon)

DupTG <- E3.TG.format %>% 
  group_by(ID_string_TG) %>% 
  summarise(n = n()) %>% 
  filter(n>1) #14 unique TG, need to remove 15

TG.clean <- E3.TG.format %>% 
  filter(ID_string_TG %in% DupTG$ID_string_TG) %>% 
  mutate(TG.avg = rowMeans(across(matches("^Y|^O")))) %>% 
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

E3.lpd.conc.fmt <- bind_rows(E3.else, E3.uniq.TG, TG.clean) #339 lipids
save(E3.lpd.conc.fmt, file = "./Output_Data/Exp3_lipid.conc.reformat.for.lipidyzer.overlap.Rdata")