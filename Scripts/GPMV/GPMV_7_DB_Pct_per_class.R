
##Double bond relative abundance per class
##Steps:
##1. Calculate class sum
##2. Take accumulative intensity of each double bond for each class divided by class sum intensity
##3. Plot double bond abundance difference between age groups only in classes that have at least 10 lipids identified
##4. Plot Phospholipid order index based on sum concentration of (PE+SM)/PC.
rm(list=ls())
library(tidyverse)
library(ggpubr)
setwd(rstudioapi::getActiveProject())

##Get class sum#############################################################################################################
load("./Output_Data/GPMV_Norm_Impt_backtoraw_all565_lipids.Rdata")
all.lipid <- raw_int.GPMV %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") %>% #4476/373=12 samples
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1))

GPMV.classsum <- all.lipid %>%
  group_by(Sample, Class) %>%
  summarise(., Class_sum = sum(Conc_Int))

save(GPMV.classsum, file = "./Output_Data/GPMV_ClassSum_565_lipids.Rdata")

##Get composition PCT for DB in each class################################################################################
load(file = "./Output_data/GPMV_DB_by_Class.Rdata")
DB_by_class_GPMV$Class <- as.character(DB_by_class_GPMV$Class)
GPMV.classsum$Class <- as.character(GPMV.classsum$Class)

GPMV_DB <- left_join(DB_by_class_GPMV, GPMV.classsum) %>%
  mutate(., DB_Pct = Sum_DB/Class_sum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(GPMV_DB, file = paste0("./Output_Data/GPMV_DB_PCT_all_samples.Rdata"))


#### DB PCT on lipids with quantitative standards####
##Get class sum#############################################################################################################
rm(list = ls())
load("./Output_Data/GPMV_Norm_Impt_log2_conc_483_lipids.Rdata")

conc.lipid <- 2^GPMV_Impt_norm_conc_all %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>% #4476/373=12 samples
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1))

GPMV.conc.classsum <- conc.lipid %>%
  group_by(Sample, Class) %>%
  summarise(., Class_sum = sum(Conc))

save(GPMV.conc.classsum, file = "./Output_Data/GPMV_ClassSum_CONC.483_lipids.Rdata")

##Get composition PCT for DB in each class################################################################################
load(file = "./Output_data/GPMV_CONC.DB_by_Class.Rdata")
CONC.DB_by_class_GPMV$Class <- as.character(CONC.DB_by_class_GPMV$Class)
GPMV.conc.classsum$Class <- as.character(GPMV.conc.classsum$Class)

GPMV_Conc.DB <- left_join(CONC.DB_by_class_GPMV, GPMV.conc.classsum, by = c("Sample", "Class")) %>%
  mutate(., DB_Pct = Sum_DB/Class_sum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(GPMV_Conc.DB, file = paste0("./Output_Data/GPMV_CONC.DB_PCT_all_samples.Rdata"))
