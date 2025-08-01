setwd(rstudioapi::getActiveProject())
library(tidyverse)
library(ggpubr)
rm(list = ls())

load("./Output_Data/Lipidyzer_qNSC_dup_TG_removed.Rdata")
Ldz.Q <- Ldz.Qui.nodup %>% 
  select(-ID_string)

Ldz.Q.org <- Ldz.Q %>%
  mutate(., ID = str_split(LipidIon, "\\.")) %>%
  rowwise() %>% 
  group_modify(~ {
    .x %>%
      mutate(., Class = ID[[1]][1]) %>%
      mutate(., Class = ifelse(grepl("DAG", Class), "DG", Class)) %>%
      mutate(., Class = ifelse(grepl("TAG", Class), "TG", Class)) %>%
      mutate(., Class = ifelse(grepl("CER", Class), "Cer", Class)) 
  }) %>% 
  relocate(Class, .after = LipidIon)

Cla.ldz.Q.nodup <- Ldz.Q.org %>% 
  rowwise() %>% 
  select(-ID) %>% 
  pivot_longer(-c(LipidIon, Class), names_to = "Samples", values_to = "Concentration") %>% 
  group_by(Samples, Class) %>% 
  summarise(., ClassSum = sum(Concentration)) %>% 
  mutate(., Age = ifelse(grepl("Y", Samples), "Young", "Old"))

save(Cla.ldz.Q.nodup, file = "./Output_Data/Lipidyzer_qNSC.ClassConcSum.dupTG.removed.Rdata")

