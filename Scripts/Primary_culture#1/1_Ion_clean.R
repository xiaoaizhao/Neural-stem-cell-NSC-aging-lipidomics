##Manually annotate Cholesterol and add it to the data matrix
##Filter lipid based on pre-determined ion for each lipid class
##444 lipids after ion clean up based on specific ion
rm(list=ls())
library(tidyverse)

setwd(rstudioapi::getActiveProject())

load("./Output_Data/0_Correct_multi_ID_by_mScore.Rdata")

##Manually label endogenous Cholesterol, labelled as "ChE(0:0)+H-H2O", with RT at around 10mins
all.id <- clean.2017.final

all.id$LipidIon[all.id$LipidIon == "ChE(0:0)+H-H2O" & grepl("10", all.id$RT_LS)] <- "Cholesterol(0:0)+H-H2O"
all.id$Class[all.id$LipidIon == "Cholesterol(0:0)+H-H2O" & grepl("10", all.id$RT_LS)] <- "Cholesterol"


#manually filled in the ion type that will be used for each lipid class, then save as csv
##import ion list to keep for each specific lipid
##Note, GM3 and LPI are two new lipid classes that were not included in prevoius ion list, manually update
ion.tbl <- read.csv("./Input_data/Final_ion_list_for_cleanUp_052620_for2017LC.csv", stringsAsFactors = F)

##add ion info to each lipid, filter to only keep one specified ion for each lipid class
Uni.clean.460 <-all.id %>%
  mutate(., Ion = paste0("[", substr(LipidIon, str_locate(LipidIon, "\\)")+1, nchar(LipidIon)), "]")) %>%
  group_by(Compound) %>%
  group_modify(~ {
    .x %>%
      mutate(., ion_match = ifelse(Ion %in% ion.tbl$Ion[ion.tbl$Class==Class], "T", "F")) %>%
      mutate(., FA_match = ifelse(str_detect(LipidIon, FattyAcid), "T", "F"))
  }) %>%
  filter(., ion_match == "T")  ##460 after ion cleanup.

Uni.clean.460 <- Uni.clean.460 %>%
  select(., -FA_match)
save(Uni.clean.460, file = paste0("./Output_data/Ion_cleaned_final_460_lipids.Rdata"))


