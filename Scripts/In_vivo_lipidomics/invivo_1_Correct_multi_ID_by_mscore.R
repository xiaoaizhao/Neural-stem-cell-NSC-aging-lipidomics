# First processing script for lipidomics on in vivo sorted NSCs

#Steps:

##1. Remove multiple IDs based on mscore, keep one compound with the highest mscore
##2. there are 30 lipids where side chain annotation in FattyAcid is different from LipidIon annotation. Therefore use annotation in FattyAcid column as side chain annotation for all lipids

#Result:
##13459 ion to go forward with ion clean in next script.
rm(list=ls())
library(tidyverse)

setwd(rstudioapi::getActiveProject())

df1 <- read.csv("./Input_Data/Xiaoai_Low_LS_th.csv", stringsAsFactors = F)

##Keep ion with unique identification
Uni.id <- df1 %>%
  filter(., !grepl("\\|", LipidIon))
multi.id <- df1 %>%
  filter(., grepl("\\|", LipidIon)) #11 lipids

##For lipids identified with multiple IDs, take one that has the highest mscore
id.df <- multi.id %>%
  mutate(., IDs = str_split(LipidIon, "\\| ")) %>%
  mutate(., scores = str_split(mscore, "\\| ")) %>%
  mutate(., RT_all = str_split(RT_LS, "\\| ")) %>%
  mutate(., FA_all = str_split(FattyAcid, "\\| ")) %>%
  group_by(Compound) %>%
  group_modify(~ {
    .x %>%
      mutate(., ID.toKeep = IDs[[1]][which.max(unlist(scores[[1]]))]) %>%
      mutate(., RT.toKeep = RT_all[[1]][which.max(unlist(scores[[1]]))]) %>%
      mutate(., FA.toKeep = FA_all[[1]][which.max(unlist(scores[[1]]))])
  }) %>%
  mutate(., LipidIon = ID.toKeep) %>%
  mutate(., RT_LS = RT.toKeep) %>%
  mutate(., FattyAcid = FA.toKeep) %>%
  select(., -matches("IDs|scores|ID.toKeep|RT_all|RT.toKeep|FA_all|FA.toKeep")) %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  ungroup() %>%
  relocate(., Compound, .after = XY.12)


####Use annotation in FattyAcid column as side chain annotation for all lipids####
exvivo.df <- bind_rows(id.df, Uni.id) %>%
  group_by(Compound) %>%
  group_modify(~ {
    .x %>%
      mutate(., LipidIon.FA = paste0(substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1),
                                     FattyAcid,
                                     substr(LipidIon, str_locate(LipidIon, "\\)")+1, nchar(LipidIon))))
  }) %>%
  arrange(., desc(LipidIon)) 

save(exvivo.df, file = paste0("./Output_Data/Multi_ID_corrected_invivo.Rdata"))




