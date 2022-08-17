##Clean duplicated lipid from 2017 LC data, 2nd processing script

##Objective: 
#1. Check to make sure all lipid has the correct number of side chains identified
#2. Remove duplicated lipids

##Steps:
#1. Make sure all lipids have the right number of side chains, result: all lipid have the correct number of sidechains - 386 total
#2. Get mean intensity of each lipid across all samples 
#3. Remove duplicated lipids (lipids that have the identical headgroup and side chain identification) by only keeping one lipid that has the highest average intensity - 373 lipids left
rm(list=ls())

library(tidyverse)
setwd(rstudioapi::getActiveProject())

load(file = "./Output_data/Ion_cleaned_final_460_lipids.Rdata")

#1 check if each lipid has the expected number of detected side chain, based on the class of lipid that it's in
##add number of side chain separation by number of "/" in ID
class_sep <- Uni.clean.460 %>%
  group_by(Class) %>%
  summarise(sep_num = max(str_count(FattyAcid, "/")))

class_sep_min <- Uni.clean.460 %>%
  group_by(Class) %>%
  summarise(sep_num = min(str_count(FattyAcid, "/")))
all.equal(class_sep,  class_sep_min) ##[1] TRUE
##Max number of side chain in each class = min number of side chain of the same class
#this is a way to show that all lipids have the correct number of detected side chain

##One more check to make sure - only keep lipid that contains the max num of side chain for each class,
##Num of lipid remains the same after filtering
sep_chain <- Uni.clean.460 %>%
  group_by(Compound) %>%
  group_modify(~ {
    .x %>%
      mutate(., complete_chain = ifelse(str_count(FattyAcid, "/") == class_sep$sep_num[class_sep$Class==Class], "T", "F"))
  }) %>%
  filter(., complete_chain == "T")  ##Still 460 after side chain cleanup.

#2 remove duplicate lipids, get mean intensity of each lipid from all samples###########
sep_chain_w_int <- sep_chain %>%
  mutate(., SideChain = str_replace(FattyAcid, "\\(", "")) %>%
  mutate(., SideChain = str_replace(SideChain, "\\)", "")) %>%
  group_modify(~ {
    .x %>%
      ungroup() %>%
      mutate(., Mean_smple_Int = mean(as.numeric(unlist(select(., matches('Y|O', ignore.case = F))))))
  })

##when encountering duplicated lipids, only keep one that has the highest average intensity from all samples, remove the rest
FA_unique373 <- sep_chain_w_int %>%
  group_by(Class, FattyAcid) %>%
  arrange(Class, FattyAcid) %>%
  group_modify(~ {
    .x %>%
      filter(., Mean_smple_Int == max(Mean_smple_Int) )
  }) ##373

##check number of unique FA for within each class and add up - total number (330) matches data frame after filtering above
check_FA_unique <- sep_chain_w_int %>%
  group_by(Class) %>%
  summarise(., nFA = length(unique(FattyAcid)))
total_num <- check_FA_unique %>%
  summarise(., num = sum(nFA)) #373, same as the number of lipids in FA_unique data.frame

save(FA_unique373, file = paste0("./Output_data/2_Ion+FA_clean_373_lipids.Rdata"))



