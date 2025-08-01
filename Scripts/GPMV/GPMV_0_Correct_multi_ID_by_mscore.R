
##First script for processing lipidomics data from GPMV
##Objective: 
#1. Process lipids with multiple identified IDs only keeping the lipid with the highest mScore.
#2. Manually annotate cholesterol
##Steps:
#1. Separate data frame into lipids that have unique ID and multiple ID
#2. Among lipids with multiple ID, only keep one identification with the highest mScore.
#3. Combine unique ID lipids with multiple ID lipids after cleaninig up
#4. Manually annotate cholesterol
rm(list=ls())
library(tidyverse)

setwd(rstudioapi::getActiveProject())

df1 <- read.csv("./Input_Data/GPMV_LS_reRun_may_2020_051820.csv", stringsAsFactors = F) ##1220 ions

##Keep ion with unique identification
Uni.id <- df1[-grep("\\|", df1$LipidIon),] ##1017 lipids
multi.id <- df1[grep("\\|", df1$LipidIon),] ##203 lipids

id.df <- multi.id %>%
  mutate(., IDs = str_split(LipidIon, "\\| ")) %>%
  mutate(., scores = str_split(mscore, "\\| ")) %>%
  mutate(., RT_all = str_split(RT_LS, "\\| ")) %>%
  group_by(Compound) %>%
  group_modify(~ {
    .x %>%
      mutate(., ID.toKeep = IDs[[1]][which.max(unlist(scores[[1]]))]) %>%
      mutate(., RT.toKeep = RT_all[[1]][which.max(unlist(scores[[1]]))])
  }) %>%
  mutate(., LipidIon = ID.toKeep) %>%
  mutate(., RT_LS = RT.toKeep) %>%
  select(., -matches("IDs|scores|ID.toKeep|RT_all|RT.toKeep")) %>%
  mutate(., FattyAcid = substr(LipidIon, str_locate(LipidIon, "\\("), str_locate(LipidIon, "\\)"))) %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1))

clean.df.GPMV <- bind_rows(id.df, Uni.id) %>%
  arrange(., as.numeric(RT_LS))

##manually annotate cholesterol
clean.df.GPMV$LipidIon[clean.df.GPMV$LipidIon == "ChE(0:0)+H-H2O" & grepl("10", clean.df.GPMV$RT_LS)] <- "Cholesterol(0:0)+H-H2O"
clean.df.GPMV$Class[grepl("Cholesterol", clean.df.GPMV$LipidIon)] <- "Cholesterol"

save(clean.df.GPMV, file = paste0("./Output_Data/GPMV_Correct_multi_ID_by_mScore.Rdata"))
