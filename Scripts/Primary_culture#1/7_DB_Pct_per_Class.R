
##Double bond relative abundance per class
##Steps:
##1. Calculate class sum
##2. Take accumulative intensity of each double bond for each class divided by class sum intensity

rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

source("./Scripts/Function_scripts/Pre-processing_functions.R")
##Get class sum####
load("./Output_Data/Spike-in_norm_Mednorm_all_373_lipids_back_to_raw_int.Rdata")
all.lipid <- raw_int %>%
  rownames_to_column(., var = "Lipid") %>%
  select(., matches("Qui|Lipid")) %>%
  pivot_longer(-Lipid, names_to = "Sample", values_to = "Intensity") %>% #4476/373=12 samples
  mutate(., Class = substr(Lipid, 1, str_locate(Lipid, "\\(")-1))

classsum <- all.lipid %>%
  group_by(Sample, Class) %>%
  summarise(., Class_sum = sum(Intensity))

save(classsum, file = "./Output_Data/Qui_NSC_ClassSum_Int_2017LC-MS.Rdata")


##Activated cell class sum####
act.all.lipid <- raw_int %>%
  rownames_to_column(., var = "Lipid") %>%
  select(., matches("Act|Lipid")) %>%
  pivot_longer(-Lipid, names_to = "Sample", values_to = "Intensity") %>% #4476/373=12 samples
  mutate(., Class = substr(Lipid, 1, str_locate(Lipid, "\\(")-1))

act.classsum <- act.all.lipid %>%
  group_by(Sample, Class) %>%
  summarise(., Class_sum = sum(Intensity))

save(act.classsum, file = "./Output_Data/Act_NSC_ClassSum_Int_2017LC-MS.Rdata")

##Get composition PCT for Qui DB in each class####
load(file = "./Output_data/Qui_DB_by_Class_2017LC.Rdata")
Qui_DB_By_Cla_2017$Class <- as.character(Qui_DB_By_Cla_2017$Class)
classsum$Class <- as.character(classsum$Class)

Qui_LC2017_df <- left_join(Qui_DB_By_Cla_2017, classsum) %>%
  mutate(., DB_Pct = Sum_DB/Class_sum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(Qui_LC2017_df, file = "./Output_Data/Qui_NSC_DB_PCT_by_Class_2017LC-MS.Rdata")

