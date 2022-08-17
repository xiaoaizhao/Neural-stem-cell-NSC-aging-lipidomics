
##Organize IS standards for Second dataset on Primary NSC culture #2

##Notes:
##Intensity from 2 different ions were collected for LPE and Cer. 
##First made sure that these two ions correlate with each other in LPE and Cer

#Steps:
##1. Remove MAG+NH4 which is not detected, Remove Cholesterol(d7) since dataset did not include cholesterol spikein
##2. Use -H for LPE and +H-H2O for Cer. Remove +H for LPE and +H for Cer
##3. Merge IS concentration with extracted intensity from each sample for calculating lipid concentration in next script
##4. Save sample IS matrix with IS concentration appended to the end.
rm(list=ls())
library(tidyverse)

setwd(rstudioapi::getActiveProject())

####Import extracted data of all IS lipids, for each sample Remove MAG+NH4 which is not detected, Remove Cholesterol(d7)####
IS <- read.csv("./Input_Data/190415_Xiaoai_SplashStandards.csv", stringsAsFactors = F, check.names = F)
IS.df <- as_tibble(IS) %>%
  select(., -contains('MAG+NH4',ignore.case = F)) %>%
  select(., -contains('Cholesterol',ignore.case = F)) %>%
  filter(., !grepl("PBS|QC", Sample)) %>%
  column_to_rownames(., var = "Sample")
IS.org <- as.data.frame(t(as.matrix(IS.df)))

####After making sure that +H vs. -H in LPE and +H vs. +H-H2O in Cer correlates very well, use -H for LPE and +H-H2O for Cer####
rem_list <- c("18_1(d7)_Lyso_PE+H","d18_1-15_0(d7)_Cer+H")
KO.IS <- IS.org %>%
  rownames_to_column(., var = "IS_class") %>%
  mutate(., Ion = paste0("[", substr(IS_class, nchar(IS_class)-1, nchar(IS_class)), "]")) %>%
  mutate(., Ion = ifelse(grepl("NH4", IS_class), "[+NH4]", Ion)) %>%
  mutate(., Ion = ifelse(grepl("+H-H2O", IS_class), "[+H-H2O]", Ion)) %>%
  filter(., !IS_class %in% rem_list) %>%
  mutate(., Class = c("LPC", "LPE", "MG", "PI", "PS", "PG", 
                      "SM", "Cer", "PC", "PE", "DG", "TG", "ChE"))

####Merge IS concentration with extracted intensity from each sample for calculating lipid concentration in next script####
labelled.conc <- read.csv("./Input_Data/KO_Labelled_IS_concentration_082720.csv", stringsAsFactors = F,
                          check.names = F, colClasses = c("character", rep("numeric", 4)))

labelled.conc.tomerge <- labelled.conc %>%
  mutate(., Class = c("PC", "LPC", "PE", "LPE", "PG", "PI", 
                      "PS", "TG", "DG", "MG", "ChE", "SM", "Cer"))
KO.IS.all <- left_join(KO.IS, labelled.conc.tomerge, by= "Class")

save(KO.IS.all, file = paste0("./Output_Data/Exp2_KO_IS_organized.Rdata"))
