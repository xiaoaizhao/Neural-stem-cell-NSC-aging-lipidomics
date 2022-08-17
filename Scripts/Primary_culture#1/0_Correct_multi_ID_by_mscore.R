

##First script for processing Primary NSC culture #1
##Objective: 
#1. Process lipids with multiple identified IDs only keeping one annotation with the highest mScore.
#2. Based on characteristic retention time for each lipid class, delete TG that elutes before 18min
library(tidyverse)

setwd(rstudioapi::getActiveProject())

df1 <- read.csv("./Input_data/2017_cell_data_Progenesis_output.csv", stringsAsFactors = F) ##843 ion

##Keep ion with unique identification
Uni.id <- df1[-grep("\\|", df1$LipidIon),] ##709 lipids
multi.id <- df1[grep("\\|", df1$LipidIon),] ##134 lipids

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

clean.df.2017 <- bind_rows(id.df, Uni.id) %>%
  arrange(., as.numeric(RT_LS))

##inspect RT range for each class
RT.summary <- clean.df.2017 %>%
  group_by(Class) %>%
  summarise(., RT_range = list(range(as.numeric(RT_LS)))) 

##Remove 3 TGs have RT <18min
clean.2017.final <- clean.df.2017 %>%
  filter(., !(Class == "TG" & as.numeric(RT_LS) < 18)) ##remove 3 TGs with RT <18 min
#840 lipids

save(clean.2017.final, file = paste0("./Output_data/0_Correct_multi_ID_by_mScore.Rdata"))
