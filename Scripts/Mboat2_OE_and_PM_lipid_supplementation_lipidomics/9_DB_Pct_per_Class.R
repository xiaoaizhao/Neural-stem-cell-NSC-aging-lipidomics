## _Mboat2_ overexpression and young plasma membrane lipid supplementation lipidomics 

##Double bond relative abundance per class
##Steps:
##1. Calculate class sum
##2. Take accumulative intensity of each double bond for each class divided by class sum intensity

rm(list = ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

##Get class sum for Mboat2 OE samples####
load("./Output_data/M2PM_backtoRAW_Norm_Impt_log2_conc_275_noconc_48_lipids.Rdata")

M2PM.all <- M2PM.raw.all

M2PM.smpl <- M2PM.all %>% #41 columns
  rownames_to_column(., var = "LipidIon") %>%
  mutate(., LipidIon = ifelse(
    grepl("Ch\\+H\\-H2O|Cholesterol", LipidIon),
    "Cholesterol(0:0)+H-H20",
    LipidIon
  )) %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Intensity") %>% #12840/321=40 samples
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1))


M2PM.classsum <- M2PM.smpl %>%
  group_by(Sample, Class) %>%
  summarise(., Class_sum = sum(Intensity))

save(M2PM.classsum, file = "./Output_Data/M2PM_all_sample_ClassSum_ConcInt.Rdata")

#=== calculate the percent concentration or intensity of DB in each class========================================================================
#=== calculate the percent concentration or intensity of DB in each class========================================================================

##Get composition PCT for Qui DB in each class####
rm(list = ls())
load(file = "./Output_Data/M2PM_all_sample_ClassSum_ConcInt.Rdata")
load(file = "./Output_data/M2PM_all_sample_DB.Rdata")


M2PM.DB.pct <- left_join(M2PM.DB_by_class, M2PM.classsum, by = c("Sample", "Class")) %>%
  mutate(., DB_Pct = Sum_DB/Class_sum) 

save(M2PM.DB.pct, file = "./Output_Data/M2PM_all_sample_DB_PCT_by_Class.Rdata")

