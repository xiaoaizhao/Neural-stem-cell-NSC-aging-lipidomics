### Double bond composition analysis

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Exp2_MsD.Norm_Impt_log2_conc_278_lipids.Rdata")

E2MsD.frmt <- 2^Exp2.log2.impt.norm.conc.MsD %>% 
  rownames_to_column(var = "LipidIon")
save(E2MsD.frmt, file = "./Output_Data/Exp2_278_conc_lipid.MSD.format.Rdata")

## ==== Aggregate lipids with the same number of double bonds for each class, in each sample ====
E2.rename <- MsD.lpd.format(E2MsD.frmt)
Exp2.DB <- E2.rename %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., Class = ifelse(LipidIon == "Cholesterol", "Cholesterol", Class)) %>% 
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) 

##Function to get aggregated intensity for each Class:DB combination####
DB_agg <- db.tally(Exp2.DB, Conc, Sample)
MsD.CONC.DB_by_class_Exp2 <- dplyr::bind_rows(DB_agg) #2496

## ==== import class sum from earlier script ====
load("./Output_Data/Exp2.KO.MsD.ClassSum.Rdata")
E2ClassSum <- Exp2.class.sum %>% 
  arrange(Class, Sample)

## ==== Calculate DB percentage across classes of each sample====
MsD.E2_CONC.DB <- left_join(MsD.CONC.DB_by_class_Exp2, E2ClassSum, by = c("Class", "Sample")) %>%
  mutate(., DB_Pct = Sum_DB/ClassSum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(MsD.E2_CONC.DB, file = "./Output_Data/MsD.Exp2.KO_CONC.DB_PCT_by_Class.Rdata")

