### Double bond composition analysis

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
source("./Scripts/Function_scripts/Pre-processing_functions.R")

load("./Output_Data/Invitro.MsD.404.log2.norm.impt.no.O8aNSC.Rdata")

Invitro.q.conc <- 2^Invitro.no.O8A.MsD %>% 
  rownames_to_column(var = "LipidIon") %>% 
  mutate(LipidIon = str_replace(LipidIon, "\\'", "")) %>% 
  select(matches("_qNSC-Q|LipidIon"))
save(Invitro.q.conc, file = "./Output_Data/Invitro_conc.lipid.qNSC.Rdata")

## ==== Aggregate lipids with the same number of double bonds for each class, in each sample ====
Invitro.rename <- MsD.lpd.format(Invitro.q.conc)
Invitro.DB <- Invitro.rename %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., Class = ifelse(LipidIon == "Cholesterol", "Cholesterol", Class)) %>% 
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) 

##Function to get aggregated intensity for each Class:DB combination####
DB_agg <- db.tally(Invitro.DB, Conc, Sample)
MsD.CONC.DB_by_class_Invitro <- dplyr::bind_rows(DB_agg) #630

## ==== import class sum from earlier script ====
load("./Output_Data//Invitro.MsD.ClassSum.A+Q.Rdata")
QClassSum <- InvitroAQ.class.sum %>% 
  filter(grepl("_qNSC-Q", Sample)) %>% 
  arrange(Class, Sample)

## ==== Calculate DB percentage across classes of each sample====
MsD.Qui_CONC.DB.invitro <- left_join(MsD.CONC.DB_by_class_Invitro, QClassSum, by = c("Class", "Sample")) %>%
  mutate(., DB_Pct = Sum_DB/ClassSum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(MsD.Qui_CONC.DB.invitro, file = "./Output_Data/Invitro.Qui_CONC.DB_PCT_by_Class.Rdata")
