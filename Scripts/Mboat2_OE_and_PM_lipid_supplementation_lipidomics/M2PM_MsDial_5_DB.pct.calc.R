### Double bond composition analysis

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
source("./Scripts/Function_scripts/Pre-processing_functions.R")

load("./Output_Data/M2PM_backtoRAW_MsD_Norm_Impt_397_lipids.Rdata") #from script #6

## ==== Aggregate lipids with the same number of double bonds for each class, in each sample ====
M2PM.n <- M2PM.MsD.lpd.conc.all %>% 
  rownames_to_column(var = "LipidIon")

M2PM.rename <- MsD.lpd.format(M2PM.n)

M2PM.DB <- M2PM.rename %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., Class = ifelse(LipidIon == "Cholesterol", "Cholesterol", Class)) %>% 
  mutate(Class = ifelse(Class == "HexCer", "Cer", Class))

### Need to fix LPC and LPE where the second side chain is listed as 0:0
LPCLPE <- M2PM.DB %>% 
  filter(Class %in% c("LPC", "LPE") & grepl("\\/", LipidIon)) %>% 
  mutate(LipidIon = paste0(substr(LipidIon, 1, str_locate(LipidIon, "\\/")-1), ")")) %>% 
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) 

othersdf <- M2PM.DB %>% 
  filter(!(Class %in% c("LPC", "LPE") & grepl("\\/", LipidIon))) %>% 
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) 

M2PM.DB.final <- bind_rows(LPCLPE, othersdf)

##Function to get aggregated intensity for each Class:DB combination####
DB_agg <- db.tally(M2PM.DB.final, Conc, Sample)
MsD.CONC.DB_by_class_M2PM <- dplyr::bind_rows(DB_agg) #2440

## ==== import class sum from earlier script ====
load("./Output_Data/M2PM.MsD.ClassSum.Rdata")

## ==== Calculate DB percentage across classes of each sample====
MsD.DB.pct.M2PM <- left_join(MsD.CONC.DB_by_class_M2PM, M2PM.class.sum, by = c("Class", "Sample")) %>%
  mutate(., DB_Pct = Sum_DB/ClassSum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(MsD.DB.pct.M2PM, file = "./Output_Data/MsD.M2PM_CONC.DB_PCT_by_Class.Rdata")
