### Double bond composition analysis

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
source("./Scripts/Function_scripts/Pre-processing_functions.R")

load("./Output_Data/Invivo.MsD.29.lipid.raw.conc.Rdata") 

## ==== Aggregate lipids with the same number of double bonds for each class, in each sample ====
Invivo.rename <- MsD.lpd.format(Invivo.raw.conc.29.good.peak)
Invivo.DB <- Invivo.rename %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., Class = ifelse(LipidIon == "Cholesterol", "Cholesterol", Class)) %>% 
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) 

##Function to get aggregated intensity for each Class:DB combination####
DB_agg <- db.tally(Invivo.DB, Conc, Sample)
MsD.CONC.DB_by_class_Invivo <- dplyr::bind_rows(DB_agg) #348

## ==== import class sum from earlier script ====
load("./Output_Data/Invivo.MsD.ClassSum.Rdata")

## ==== Calculate DB percentage across classes of each sample====
MsD.DB.pct.invivo <- left_join(MsD.CONC.DB_by_class_Invivo, Invivo.class.sum, by = c("Class", "Sample")) %>%
  mutate(., DB_Pct = Sum_DB/ClassSum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(MsD.DB.pct.invivo, file = "./Output_Data/MsD.Invivo_CONC.DB_PCT_by_Class.Rdata")

