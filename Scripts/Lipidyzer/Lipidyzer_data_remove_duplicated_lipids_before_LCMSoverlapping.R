### Remove duplicated TG lipids in Lipidyzer data before comparing with other datasets

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(eulerr)
source("./Scripts/Function_scripts/Effect_size_functions.R")

load("./Output_data/Ldz_511_lipids_backtoraw.conc_MedConc_norm_imputed.Rdata")
Dyzer.raw <- ldz.raw_conc %>% 
  rownames_to_column(var = "LipidIon") %>% 
  mutate(., Class = substr(LipidIon, 1, stri_locate_first(LipidIon, regex =  "\\.")-1)) %>% 
  mutate(., Class = case_when(
    grepl("^TAG", LipidIon) ~ "TG",
    grepl("^DAG", LipidIon) ~ "DG",
    !grepl("^TAG | ^DAG", LipidIon) ~ Class
  ))

dyzer.cla <- unique(Dyzer.raw$Class)

## ====Remove lipid classes that were not detected in LC-MS====
# Remove classes where data can't be compared to LC-MS
# FFA is not quantified, SM is only reported with one chain in lipidyzer and 
# Cer is only reported with one chain
# Re-format side chains

Dyzer.ClaFltr <- Dyzer.raw %>% #447 lipids
  filter(!Class %in% c("HCER", "LCER", "DCER", "FFA", "CER", "SM")) %>% 
  rowwise() %>% 
  mutate(., n.str = c(str_split(LipidIon, "\\."))) %>% 
  mutate(LipidIon1 = LipidIon) %>% 
  group_by(LipidIon1) %>% 
  group_modify(~{
    .x %>% 
  mutate(., ID_string = case_when(
    Class %in% c("CE", "LPC", "LPE") ~ list(c(n.str[[1]][1], paste0(n.str[[1]][2],":", n.str[[1]][3]))),
    Class %in% c("DG", "PC", "PE") ~ list(c(n.str[[1]][1], paste0(n.str[[1]][2],":", n.str[[1]][3]), paste0(n.str[[1]][4],":", n.str[[1]][5]))),
    Class %in% c("TG") ~ list(c("TG", paste0(substr(n.str[[1]][1], 4, 5),":", n.str[[1]][2])))))
    })%>% 
  mutate(., ID_string = ifelse(Class == "DG", list(c("DG", ID_string[[1]][2], ID_string[[1]][3])), ID_string)) %>% 
  mutate(., ID_string = ifelse(grepl("PE.P|PE.O", LipidIon), 
                               list(c(n.str[[1]][1], paste0(n.str[[1]][2],"-", n.str[[1]][3], ":", n.str[[1]][4]),
                                      paste0(n.str[[1]][5],":", n.str[[1]][6]))),
                               ID_string)) %>% 
  mutate(., ID_string = ifelse(Class == "CE", list(c("ChE", paste0(n.str[[1]][2],":", n.str[[1]][3]))), ID_string)) %>% 
  relocate(c(n.str, ID_string), .after = "LipidIon") 

## ====Quiescent cells====
### Qui remove duplicated TG from lipidyzer data
##After - removed 150 duplicated TG

Qui.Dyzer <- Dyzer.ClaFltr %>% 
  ungroup() %>% 
  select( "LipidIon"|"ID_string"|matches("._qui")) 

qui.mean <- Qui.Dyzer %>%
  select(-ID_string) %>% 
  column_to_rownames(var = "LipidIon") %>% 
  mutate(Lpd.Mean = rowMeans(.)) %>% 
  rownames_to_column(var = "LipidIon")

Qui.nodup.ls <- left_join(Qui.Dyzer, qui.mean) %>% 
  mutate(ID_string1 = ID_string) %>% 
  group_by(ID_string1) %>% 
  group_modify(~{
    .x %>% 
    filter(Lpd.Mean == max(Lpd.Mean))
  }) %>% 
  ungroup() %>% 
  select(ID_string, Lpd.Mean, LipidIon)


#### Save qNSC samples from Lipidzyer data with this subset of lipids

Ldz.Qui.nodup <- Qui.Dyzer %>% 
  filter(LipidIon %in% Qui.nodup.ls$LipidIon)

save(Ldz.Qui.nodup, file = "./Output_Data/Lipidyzer_qNSC_for.LC-MS.ovlp.Rdata")

## ====Activated cells====
Act.Dyzer <- Dyzer.ClaFltr %>% 
  ungroup() %>% 
  select( "LipidIon"|"ID_string"|matches("._act")) 

Act.mean <- Act.Dyzer %>%
  select(-ID_string) %>% 
  column_to_rownames(var = "LipidIon") %>% 
  mutate(Lpd.Mean = rowMeans(.)) %>% 
  rownames_to_column(var = "LipidIon")

Act.nodup.ls <- left_join(Act.Dyzer, Act.mean) %>% 
  mutate(ID_string1 = ID_string) %>% 
  group_by(ID_string1) %>% 
  group_modify(~{
    .x %>% 
      filter(Lpd.Mean == max(Lpd.Mean))
  }) %>% 
  ungroup() %>% 
  select(ID_string, Lpd.Mean, LipidIon)


#### Save qNSC samples from Lipidzyer data with this subset of lipids

Ldz.Act.nodup <- Act.Dyzer %>% 
  filter(LipidIon %in% Act.nodup.ls$LipidIon)

save(Ldz.Act.nodup, file = "./Output_Data/Lipidyzer_aNSC_for.LC-MS.ovlp.Rdata")