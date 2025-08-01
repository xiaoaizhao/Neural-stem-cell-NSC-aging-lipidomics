
##Double bond relative abundance per class
##Steps:
##1. Calculate class sum
##2. Take accumulative intensity of each double bond for each class divided by class sum intensity
##The above 2 steps are identical to Primary culture LC-MS data
##3. Plot double bond abundance difference between age groups in all classes rather than filter by class size - 
##This is because in vivo data has significantly less number of detected lipids. With a class size filter >10 there will only be 4 classes left.
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Pre-processing_functions.R")

####==== All lipids ====
##Get class sum
load("./Output_Data/Invivo_Norm_Impt_backtoraw_all130_lipids.Rdata")
all.lipid <- Invivo.raw_int.invivo %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") %>% #4476/373=12 samples
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1))

Invivo.classsum <- all.lipid %>%
  group_by(Sample, Class) %>%
  summarise(., Class_sum = sum(Conc_Int))

save(Invivo.classsum, file = "./Output_Data/Invivo_Class.sum_ConcInt.Rdata")
##Get composition PCT for DB in each class
load(file = "./Output_data/Invivo_DB_by_Class_qNSC.Rdata")
DB_by_class_invivo$Class <- as.character(DB_by_class_invivo$Class)
Invivo.classsum$Class <- as.character(Invivo.classsum$Class)

Invivo.DB.PCT <- left_join(DB_by_class_invivo, Invivo.classsum) %>%
  mutate(., DB_Pct = Sum_DB/Class_sum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(Invivo.DB.PCT, file = "./Output_Data/Invivo_DB_PCT_by_Class.Rdata")

####==== Lipids with quantitative standards  ====
##Get class sum
load("./Output_Data/Invivo_Norm_Impt_log2_conc_121_lipids.Rdata")
conc.lipid <- 2^Invivo.Impt_norm_conc_all %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>% #4476/373=12 samples
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1))

Invivo.CONC.classsum <- conc.lipid %>%
  group_by(Sample, Class) %>%
  summarise(., Class_sum = sum(Conc))

Invivo.CONC.lipidN <- conc.lipid %>% 
  group_by(Sample, Class) %>% 
  tally() %>% 
  group_by(Class) %>% 
  summarise(., AvgNlpd.per.class = mean(n))

save(Invivo.CONC.classsum, file = "./Output_Data/Invivo_CONC.Class.sum.Rdata")

##Get composition PCT for DB in each class####
load(file = "./Output_data/Invivo_CONC.DB_by_Class_qNSC.Rdata")
DB_CONC.by_class_invivo$Class <- as.character(DB_CONC.by_class_invivo$Class)
Invivo.CONC.classsum$Class <- as.character(Invivo.CONC.classsum$Class)

Invivo.CONC.DB.PCT <- left_join(DB_CONC.by_class_invivo, Invivo.CONC.classsum) %>%
  mutate(., DB_Pct = Sum_DB/Class_sum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(Invivo.CONC.DB.PCT, file = "./Output_Data/Invivo_CONC_DB_PCT_by_Class.Rdata")


