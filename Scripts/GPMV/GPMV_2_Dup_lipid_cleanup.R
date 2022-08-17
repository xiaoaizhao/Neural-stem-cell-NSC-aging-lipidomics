##Clean duplicated lipid from GPMV data, 2nd processing script

##Objective: 
#1. Check to make sure all lipid has the correct number of side chains identified
#2. Remove duplicated lipids

##Steps:
#1. Make sure all lipids have the right number of side chains, result: all lipid have the correct number of sidechains - 728 total
#2. Get mean intensity of each lipid across all samples 
#3. Remove duplicated lipids (lipids that have the identical headgroup and side chain identification) by only keeping one lipid that has the highest average intensity - 565 lipids left
rm(list=ls())
library(tidyverse)

setwd(rstudioapi::getActiveProject())

load(file = "./Output_Data/GPMV_Ion_cleaned_final_728_lipids.Rdata")

#1 check if each lipid in has the expected number of detected side chain, based on the class of lipid that it's in####
##add number of side chain separation by number of "/" in ID
class_sep <- GPMV.Uni.clean.728 %>%
  group_by(Class) %>%
  summarise(sep_num = max(str_count(FattyAcid, "/")))

class_sep_min <- GPMV.Uni.clean.728 %>%
  group_by(Class) %>%
  summarise(sep_num = min(str_count(FattyAcid, "/")))
all.equal(class_sep,  class_sep_min) ##[1] TRUE
##Max number of side chain in each class = min number of side chain of the same class
#this is a way to show that all lipids have the correct number of detected side chain

##One more check to make sure - only keep lipid that contains the max num of side chain for each class,
##Num of lipid remains the same after filtering
sep_chain <- GPMV.Uni.clean.728 %>%
  group_by(Compound) %>%
  group_modify(~ {
    .x %>%
      mutate(., complete_chain = ifelse(str_count(FattyAcid, "/") == class_sep$sep_num[class_sep$Class==Class], "T", "F"))
  }) %>%
  filter(., complete_chain == "T")  ##Still 728 after side chain cleanup.

#2 remove duplicate lipids, get mean intensity of each lipid from all samples####
t <- sep_chain %>%
  mutate(., SideChain = str_replace(FattyAcid, "\\(", "")) %>%
  mutate(., SideChain = str_replace(SideChain, "\\)", "")) %>%
  group_modify(~ {
    .x %>%
      ungroup() %>%
      mutate(., Mean_smple_Int = mean(as.numeric(unlist(select(., matches('Y|O', ignore.case = F))))))
  })

##when encountering duplicated lipids, only keep one that has the highest average intensity from all samples, remove the rest
GPMV.FA_unique565 <- t %>%
  group_by(Class, FattyAcid) %>%
  arrange(Class, FattyAcid) %>%
  group_modify(~ {
    .x %>%
      filter(., Mean_smple_Int == max(Mean_smple_Int) )
  }) ##565

check_FA_unique <- t %>%
  group_by(Class) %>%
  summarise(., nFA = length(unique(FattyAcid)))

total_num <- check_FA_unique %>%
  summarise(., num = sum(nFA)) #565, same as the number of lipids in FA_unique data.frame
total_num

save(GPMV.FA_unique565, file = paste0("./Output_Data/GPMV_Ion+FA_clean_565_lipids.Rdata"))

