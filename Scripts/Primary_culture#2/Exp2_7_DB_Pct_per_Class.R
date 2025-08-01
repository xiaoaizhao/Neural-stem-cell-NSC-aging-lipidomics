
##Double bond relative abundance per class
##Steps:
##1. Calculate class sum
##2. Take accumulative intensity of each double bond for each class divided by class sum intensity
##3. Plot double bond abundance difference between age groups only in classes that have at least 10 lipids identified
rm(list=ls())
library(tidyverse)
source("./Scripts/Function_scripts/Effect_size_functions.R")
setwd(rstudioapi::getActiveProject())

##Get class sum#############################################################################################################
load("./Output_Data/Exp2_Norm_Impt_backtoraw_all693_lipids.Rdata")
all.lipid <- raw_int.exp2 %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Intensity") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., KO = substr(Sample, str_locate(Sample, "_"), nchar(Sample))) %>%
  mutate(., ID = substr(Sample, 1, str_locate(Sample, "_")-1))

Exp2.classsum <- all.lipid %>%
  group_by(KO, ID, Class) %>%
  summarise(., Class_sum = sum(Conc_Intensity))
save(Exp2.classsum, file = paste0("./Output_Data/Exp2_class_sum_693_lipids.Rdata"))

##Get composition PCT for DB in each class################################################################################
load("./Output_Data/Exp2_DB_by_Class.Rdata")
Exp2_DB_by_class$Class <- as.character(Exp2_DB_by_class$Class)
Exp2.classsum$Class <- as.character(Exp2.classsum$Class)

Exp2.classsum <- Exp2.classsum %>%
  mutate(., Sample = paste0(ID, KO))

Exp2_DB <- left_join(Exp2_DB_by_class, Exp2.classsum) %>%
  mutate(., DB_Pct = Sum_DB/Class_sum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(Exp2_DB, file = paste0("./Output_Data/Exp2_DB_PCT_all_samples.Rdata"))

#### DB PCT analysis on lipids with quantitative standards only ####
##Get class sum#############################################################################################################
load("./Output_Data/Exp2_Norm_Impt_backtoraw_conc612_lipids.Rdata")
conc.lipid <- raw_conc.exp2 %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., KO = substr(Sample, str_locate(Sample, "_"), nchar(Sample))) %>%
  mutate(., ID = substr(Sample, 1, str_locate(Sample, "_")-1))

Exp2.CONC.classsum <- conc.lipid %>%
  group_by(KO, ID, Class) %>%
  summarise(., Class_sum = sum(Conc))
save(Exp2.CONC.classsum, file = paste0("./Output_Data/Exp2_class_sum_CONC.612_lipids.Rdata"))

##Get composition PCT for DB in each class################################################################################
load("./Output_Data/Exp2_CONC.DB_by_Class.Rdata")
Exp2_CONC.DB_by_class$Class <- as.character(Exp2_CONC.DB_by_class$Class)
Exp2.CONC.classsum$Class <- as.character(Exp2.CONC.classsum$Class)

Exp2.CONC.classsum <- Exp2.CONC.classsum %>%
  mutate(., Sample = paste0(ID, KO))

Exp2_CONC.DB <- left_join(Exp2_CONC.DB_by_class, Exp2.CONC.classsum) %>%
  mutate(., DB_Pct = Sum_DB/Class_sum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(Exp2_CONC.DB, file = paste0("./Output_Data/Exp2_CONC.DB_PCT_all_samples.Rdata"))