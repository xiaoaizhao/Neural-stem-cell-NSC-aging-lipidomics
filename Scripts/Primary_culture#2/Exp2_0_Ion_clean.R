
##Primary NSC culture #2, this set of samples were transduced with gDNA either targeting specific genes or non-targeting as control

####Steps:
##1. Make sure all lipids are identified with unique ID
##2. Manually label endogenous cholesterol and add it to the data matrix
##3. Filter lipid based on pre-determined ion for each lipid class, exactly as done in the 1st exp on primary NSC culture by LC-MS
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

df <- read.csv("./Input_Data/190411_Xiaoai_lipids_final_w_score.csv", stringsAsFactors = F)

Uni.id <- df[!grepl("\\|", df$LipidIon),] ##2143 lipids, every single lipid is uniquely identified

##Manually inspect all ChE ion, Cholesterol is ChE(0:0)+H-H2O that eludes at around 10mins
ChE <- Uni.id %>%
  filter(., grepl("ChE", LipidIon))

Uni.id$LipidIon[Uni.id$LipidIon == "ChE(0:0)+H-H2O" & grepl("9", Uni.id$RT)] <- "Cholesterol(0:0)+H-H2O"
Uni.id$Class[grepl("Cholesterol", Uni.id$LipidIon)] <- "Cholesterol"

##import ion list to keep for each specific lipid, this was determined from previous experiment based on the most abundant ion of each lipid class
ion.tbl <- read.csv("./Input_Data/Final_ion_list_for_cleanUp_052620_for2017LC.csv", stringsAsFactors = F)


##only include 28 classes that have pre-determined ion for quantification
KO.clean.class.filter <- Uni.id %>%
  filter(., Class %in% ion.tbl$Class) #2114

##add ion info to each lipid
KO.clean<-KO.clean.class.filter %>%
  mutate(., Ion = paste0("[", substr(LipidIon, str_locate(LipidIon, "\\)")+1, nchar(LipidIon)), "]")) %>%
  mutate(., Compound = paste(LipidIon, RT, sep = "_")) %>%
  group_by(., Compound) %>%
  group_modify(~ {
    .x %>%
      mutate(., ion_match = ifelse(Ion == ion.tbl$Ion[ion.tbl$Class==Class], "T", "F"))
  }) %>%
  filter(., ion_match == "T") %>%  
  arrange(., Ion) ##1114 after ion cleanup.

save(KO.clean, file = paste0("./Output_Data/Exp2_Ion_clean_1114_lipids.Rdata"))
