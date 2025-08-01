
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
source("./Scripts/Function_scripts/Pre-processing_functions.R")

load("./Output_Data/Exp3_qNSC_quant.lipids.Rdata")
# ==== Lipids with quantitative standards====
## ==== Aggregate lipids with the same number of double bonds for each class, in each sample ====
Exp3.DB <- Exp3.conc.Q %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) 

##Function to get aggregated intensity for each Class:DB combination####
DB_agg <- db.tally(Exp3.DB, Conc, Sample)
CONC.DB_by_class_Exp3 <- dplyr::bind_rows(DB_agg) #1608

save(CONC.DB_by_class_Exp3, file = paste0("./Output_Data/Exp3_Qui.CONC.DB.Rdata"))

## ==== Calculate class sum concentration of each sample====
Exp3.class.sum <- Exp3.DB %>% 
  group_by(Class, Sample) %>% 
  summarise(ClassSum = sum(Conc))

## ==== Calculate DB percentage across classes of each sample====
Qui_CONC.DB.E3 <- left_join(CONC.DB_by_class_Exp3, Exp3.class.sum, by = c("Class", "Sample")) %>%
  mutate(., DB_Pct = Sum_DB/ClassSum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(Qui_CONC.DB.E3, file = "./Output_Data/Exp3.Qui_CONC.DB_PCT_by_Class.Rdata")

# ==== All lipids ====
load("./Output_Data/Exp3_all.lipid.all.cells.Rdata")

## ==== Aggregate lipids with the same number of double bonds for each class, in each sample ====
Exp3.DB.all <- Exp3.all.lpd.all.cell %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) 

##Function to get aggregated intensity for each Class:DB combination####
DB_agg <- db.tally(Exp3.DB.all, Conc_Int, Sample)
DB_by_class_Exp3 <- dplyr::bind_rows(DB_agg) #1608
  
save(DB_by_class_Exp3, file = paste0("./Output_Data/Exp3_all.DB.all.sample.Rdata"))

## ==== Calculate class sum concentration of each sample====
Exp3.class.sum <- Exp3.DB.all %>% 
  group_by(Class, Sample) %>% 
  summarise(ClassSum = sum(Conc_Int))

## ==== Calculate DB percentage across classes of each sample====
ALL.DB.E3 <- left_join(DB_by_class_Exp3, Exp3.class.sum, by = c("Class", "Sample")) %>%
  mutate(., DB_Pct = Sum_DB/ClassSum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

Qui.DB.E3 <- ALL.DB.E3 %>% 
  filter(grepl("_qNSC-Q", Sample))
save(Qui.DB.E3, file = "./Output_Data/Exp3.Qui_all.DB_PCT_by_Class.Rdata")
