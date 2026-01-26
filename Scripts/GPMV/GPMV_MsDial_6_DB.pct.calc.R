### Double bond composition analysis

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
source("./Scripts/Function_scripts/Pre-processing_functions.R")

load("./Output_Data/GPMV_MsD.Raw.conc.Norm_Impt_291_lipids.Rdata") #from script #3

## ==== Aggregate lipids with the same number of double bonds for each class, in each sample ====
GPMV.n <- GPMV.impt.norm.raw.conc.lpd.MsD %>% 
  rownames_to_column(var = "LipidIon")

GPMV.rename <- MsD.lpd.format(GPMV.n)

GPMV.DB <- GPMV.rename %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., Class = ifelse(LipidIon == "Cholesterol", "Cholesterol", Class)) %>% 
  mutate(Class = ifelse(Class == "HexCer", "Cer", Class))

### Need to fix LPC and LPE where the second side chain is listed as 0:0
LPCLPE <- GPMV.DB %>% 
  filter(Class %in% c("LPC", "LPE") & grepl("\\/", LipidIon)) %>% 
  mutate(LipidIon = paste0(substr(LipidIon, 1, str_locate(LipidIon, "\\/")-1), ")")) %>% 
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) 

othersdf <- GPMV.DB %>% 
  filter(!(Class %in% c("LPC", "LPE") & grepl("\\/", LipidIon))) %>% 
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) 

GPMV.DB.final <- bind_rows(LPCLPE, othersdf)

##Function to get aggregated intensity for each Class:DB combination####
DB_agg <- db.tally(GPMV.DB.final, Conc, Sample)
MsD.CONC.DB_by_class_GPMV <- dplyr::bind_rows(DB_agg) #896

## ==== import class sum from earlier script ====
load("./Output_Data/GPMV.MsD.ClassSum.Rdata")

## ==== Calculate DB percentage across classes of each sample====
MsD.DB.pct.GPMV <- left_join(MsD.CONC.DB_by_class_GPMV, GPMV.class.sum, by = c("Class", "Sample")) %>%
  mutate(., DB_Pct = Sum_DB/ClassSum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(MsD.DB.pct.GPMV, file = "./Output_Data/MsD.GPMV_CONC.DB_PCT_by_Class.Rdata")
