
##4th Script on in vivo sorted NSC lipidomics dataset

##Steps:
#1. Make sure all lipids have the right number of side chains, result: 14 phospholipids only have one side chain reported
#2. After inspecting the spectra, manually recovered 3 phospholipids, annotated side chain info - 148 lipids total
#3. Remove in each class that have the exact same FA (regardless of total carbon/double bond number of sn position) - 130 lipids left. Delete 18 duplicated lipids
rm(list=ls())
library(tidyverse)

setwd(rstudioapi::getActiveProject())

load(file = "./Output_Data/Rename_sample_ion_clean_159_lipids.Rdata")
#1 check if each lipid in has the expected number of detected side chain, based on the class of lipid that it's in
##add number of side chain separation by number of "/" in ID
class_sep <- rename.clean %>%
  group_by(Class) %>%
  summarise(sep_num = max(str_count(FattyAcid, "_")))

sep_chain <- rename.clean %>%
  group_by(Compound) %>%
  group_modify(~ {
    .x %>%
      mutate(., complete_chain = ifelse(str_count(FattyAcid, "_") == class_sep$sep_num[class_sep$Class==Class], "T", "F"))
  }) %>%
  filter(., complete_chain == "T")  ##14 phospholipid only have one side chain reported

##Looks like there are PC and PE and PS that don't have side chain separated, after inspecting their spectra, recovered 3 phospholipids####
check.list <- c("PC", "PE", "PS")
PCPEPS_noSep <- rename.clean %>%
  filter(., Class %in% check.list) %>%
  filter(., !grepl("_", LipidIon.FA))


#keep.list <- c("9.99_732.5524m/z", "9.85_808.5848m/z", "11.71_810.5272m/z")

Keep1 <- PCPEPS_noSep %>%
  filter(., Compound == "9.99_732.5524m/z") %>%
  mutate(., LipidIon.FA = "PC(16:0_16:1)")
Keep2 <- PCPEPS_noSep %>%
  filter(., Compound == "9.85_808.5848m/z") %>%
  mutate(., LipidIon.FA = "PC(16:0_22:5)")
Keep3 <- PCPEPS_noSep %>%
  filter(., Compound == "11.71_810.5272m/z") %>%
  mutate(., LipidIon.FA = "PS(18:1_20:3)")

Phospho_toKeep <- bind_rows(Keep1, Keep2, Keep3) %>%
  mutate(., FattyAcid = substr(LipidIon.FA, str_locate(LipidIon.FA, "\\("), str_locate(LipidIon.FA, "\\)")))

Sep_chain_all <- bind_rows(sep_chain %>% select(., -complete_chain), Phospho_toKeep)  #148 lipids with correct side chain

#2 remove duplicate lipids, get mean intensity of each lipid from all samples
t <- Sep_chain_all %>%
  mutate(., SideChain = str_replace(FattyAcid, "\\(", "")) %>%
  mutate(., SideChain = str_replace(SideChain, "\\)", "")) %>%
  group_modify(~ {
    .x %>%
      ungroup() %>%
      mutate(., Mean_smple_Int = mean(as.numeric(unlist(select(., matches('Y|O', ignore.case = F))))))
  })

FA_unique130 <- t %>%
  group_by(Class, FattyAcid) %>%
  arrange(Class, FattyAcid) %>%
  group_modify(~ {
    .x %>%
      filter(., Mean_smple_Int == max(Mean_smple_Int) )
  }) ##130

check_FA_unique <- t %>%
  group_by(Class) %>%
  summarise(., nFA = length(unique(FattyAcid)))

total_num <- check_FA_unique %>%
  summarise(., num = sum(nFA)) #130, same as the number of lipids in FA_unique data.frame
total_num


save(FA_unique130, file = paste0("./Output_Data/Invivo_Ion+FA_clean_130_lipids.Rdata"))

