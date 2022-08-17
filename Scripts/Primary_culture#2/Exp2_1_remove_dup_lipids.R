
##Clean duplicated lipid from Second dataset on Primary NSC culture #2
#1 make sure all lipids have the right number of side chains, result: all lipid have the correct number of sidechains - 1114 total
#2 remove in each class that have the exact same FA (regardless of total carbon/double bond number of sn position) - 738 lipids left
#Note on step 2: Upon inspection, found that there are 2 PC species with the identical side chain annotation AND intensity in all samples
#Remove one to only keep one for each species with a slightly later RT.
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

load(file = "./Output_Data/Exp2_Ion_clean_1114_lipids.Rdata")
#1 check if each lipid has the expected number of detected side chain, based on the class of lipid that it's in
##add number of side chain separation by number of "/" in ID
class_sep <- KO.clean %>%
  group_by(Class) %>%
  summarise(sep_num = max(str_count(FattyAcid, "/")))

class_sep_min <- KO.clean %>%
  group_by(Class) %>%
  summarise(sep_num = min(str_count(FattyAcid, "/")))
all.equal(class_sep,  class_sep_min) ##[1] TRUE, this is another way to show that all lipids have the correct chain length

#2 remove duplicate lipids, get mean intensity of each lipid from all samples###########
t <- KO.clean %>%
  mutate(., SideChain = str_replace(FattyAcid, "\\(", "")) %>%
  mutate(., SideChain = str_replace(SideChain, "\\)", "")) %>%
  group_modify(~ {
    .x %>%
      ungroup() %>%
      mutate(., Mean_smple_Int = mean(as.numeric(unlist(select(., matches('Y|O', ignore.case = F))))))
  })

##when encountering duplicated lipids, only keep one that has the highest average intensity from all samples, remove the rest
KO.FA_unique <- t %>%
  group_by(Class, FattyAcid) %>%
  arrange(Class, FattyAcid) %>%
  group_modify(~ {
    .x %>%
      filter(., Mean_smple_Int == max(Mean_smple_Int) )
  }) ##740

##Upon inspection, found that there are 2 PC species with the identical side chain annotation AND intensity in all samples
##Remove one to only keep one for each species with a slightly later RT.
PC <- KO.FA_unique %>%
  filter(., Class == "PC") #198
dupPC <- PC[grep("17:1/22:6|15:0/16:0",PC$FattyAcid),]
dupPC$Mean_smple_Int #[1] 1742500000 1742500000   34625000   34625000

remPC <- dupPC$Compound[c(1,3)] ##these two have slightly earlier RT, will remove from final matrix to keep the other two that have slightly later RT

KO.FA_clean_738 <- KO.FA_unique %>% ##remove 2 dup PC
  filter(., !Compound %in% remPC)

save(KO.FA_clean_738, file = paste0("./Output_Data/Exp2_FA_clean_738_lipids.Rdata"))
