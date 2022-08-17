
##Double bond relative abundance per class
##Steps:
##1. Calculate class sum
##2. Take accumulative intensity of each double bond for each class divided by class sum intensity
##3. Plot double bond abundance difference between age groups only in classes that have at least 10 lipids identified
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

##Get composition PCT for Act DB in each class####
load("./Output_Data/Act_NSC_ClassSum_Int_2017LC-MS.Rdata")
load("./Output_Data/Act_Qui_DB_By_Class_2017LC.Rdata")
Act_DB <- DB_by_class_2017 %>%
  filter(., grepl("Act", Sample, ignore.case = F))

Act_LC2017_df <- left_join(Act_DB, act.classsum) %>%
  mutate(., DB_Pct = Sum_DB/Class_sum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(Act_LC2017_df, file = "./Output_Data/Act_NSC_DB_PCT_by_Class_2017LC-MS.Rdata")

##DB pct age comparison, fold change between old and young####
OY_FC <- DB.FC(Qui_LC2017_df)

####Filter to only plot class that contains at least 10 lipids####
Lipid_numLC <- raw_int %>%
  mutate(., LipidIon = rownames(raw_int)) %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  group_by(Class) %>%
  summarise(., nLipid = n()) %>%
  filter(.,nLipid>=10 )

LC_DB_filter <- OY_FC %>%
  filter(., Class %in% Lipid_numLC$Class)

##plot with class list which has at least 10 lipids detected in this experiments#########################################
##plot with class list which has at least 10 lipids detected in this experiments#########################################
p <- ggplot(data =  LC_DB_filter, aes(x= Class, y= DB_num, fill= log2FC_OvY)) 
p +  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0)+
  labs(fill= "Log2 FC")+
  labs(title = "2017 LC Qui Old vs. Young DB abundance", x = "Class", y = "Double Bond Number", color = "")+
  theme(text=element_text(size = 13, face = "plain"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(paste0("./Figure_Panels/Fig_1e.pdf"), width = 5, height = 5, useDingbats=FALSE)

