
##overlap LC-MS with Lipidyzer
##from common detected class of lipids clas_commmon_w_lc <- c("DG", "LPC", "LPE", "PC", "PE", "TG")
##Following script "Organize_Lipidyzer_LC-MS_for_Overlap.R"
##remove duplicates in both Lipidyzer and LC-MS dataset
##Find common lipids in both datasets
##Plot with UpSet plot to show common +unique lipids from both datasets
rm(list=ls())
library(tidyverse)

setwd(rstudioapi::getActiveProject())
load("./Output_data/LC-MS_FA1_for_dyzer_overlap.Rdata")
load("./Output_data/LC-MS_FA2_for_dyzer_overlap.Rdata")
load("./Output_data/Lipidyzer_FA1_for_overlap.Rdata")
load("./Output_data/Lipidyzer_FA2_for_overlap.Rdata")

####remove duplicates in both Lipidyzer and LC-MS dataset####
Dyzer <- rbind(Dyzer.FA1, Dyzer.FA2)

Group.list <- Dyzer %>%
  ungroup() %>%
  group_by(ID_string) %>%
  summarise(., group_size = n()) 

##list of duplicated TG
dup.list <- Group.list$ID_string[Group.list$group_size>1] #43
single.list <- Group.list$ID_string[Group.list$group_size==1] #249

##Remove duplicates by keeping the one that has the highest effect size (absolute value)
Dyzer.rmv.dup <- Dyzer %>%
  ungroup() %>%
  filter(., ID_string %in% dup.list) %>%
  group_by(ID_string) %>%
  group_modify(~ {
    .x %>%
      filter(., abs(es_g) == max(abs(es_g)))
  }) 

Dyzer.single <- Dyzer %>%
  filter(., ID_string %in% single.list)
Dyzer.no.dup <- bind_rows(Dyzer.rmv.dup, Dyzer.single) #292 442-292=150 duplicated TG

save(Dyzer.no.dup, file = paste0("./Output_data/Lipidyzer_remove_TG_dup.Rdata"))

####LC-MS
LC <- rbind(LC.FA1.all, LC.FA2.all)
LC.no.dup<- LC %>%
  group_by(ID_string) %>%
  group_modify(~ {
    .x %>%
      filter(., abs(es_g) == max(abs(es_g)))
  }) #229, instead of 232, removed 3 duplicated lipids
save(LC.no.dup, file = paste0("./Output_data/LC_remove_dup.Rdata"))


############################Overlap############################################################
load("./Output_data/LC_remove_dup.Rdata")
load("./Output_data/Lipidyzer_remove_TG_dup.Rdata")
LC.in.Dyzer <- LC.no.dup

Dyzer.in.LC <- Dyzer.no.dup

Dyzer.slim <- Dyzer.in.LC %>%
  select(., matches("ID_string|Lipidyzer|LipidIon")) #292

LC.slim <- LC.in.Dyzer %>%
  ungroup() %>%
  select(., matches("ID_string|LC|LipidIon", ignore.case = F)) #229

LC_dyzer_overlap <- inner_join(LC.slim, Dyzer.slim, by = "ID_string")
save(LC_dyzer_overlap, file = paste0("./Output_data/LC_Lipidyzer_overlap_115_lipids_from_6_clas.Rdata" ))
