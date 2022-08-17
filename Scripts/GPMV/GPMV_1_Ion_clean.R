
##Filter lipid based on pre-determined ion for each lipid class
##728 lipids after ion clean up based on specific ion
rm(list=ls())
library(tidyverse)

setwd(rstudioapi::getActiveProject())
load("./Output_Data/GPMV_Correct_multi_ID_by_mScore.Rdata")

##import ion list to keep for each specific lipid, this was determined from previous experiment based on the most abundant ion of each lipid class
ion.tbl <- read.csv("./Input_Data/Final_ion_list_for_cleanUp_052620_for2017LC.csv", stringsAsFactors = F)

##add ion info to each lipid
Uni.clean.728 <-clean.df.GPMV %>%
  mutate(., Ion = paste0("[", substr(LipidIon, str_locate(LipidIon, "\\)")+1, nchar(LipidIon)), "]")) %>%
  group_by(., Compound) %>%
  group_modify(~ {
    .x %>%
      mutate(., ion_match = ifelse(Ion == ion.tbl$Ion[ion.tbl$Class==Class], "T", "F")) %>%
      mutate(., FA_match = ifelse(str_detect(LipidIon, coll(FattyAcid)), "T", "F"))
  }) %>%
  filter(., ion_match == "T") %>%  ##728 after ion cleanup.
  arrange(., Ion)

GPMV.Uni.clean.728 <- Uni.clean.728 %>%
  select(., -matches("FA_match|ion_match"))

save(GPMV.Uni.clean.728, file = paste0("./Output_Data/GPMV_Ion_cleaned_final_728_lipids.Rdata"))
