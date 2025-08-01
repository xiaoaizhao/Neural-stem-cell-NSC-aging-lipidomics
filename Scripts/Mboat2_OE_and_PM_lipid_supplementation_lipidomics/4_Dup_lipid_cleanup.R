##Clean duplicated lipid from _Mboat2_ overexpression and young plasma membrane lipid supplementation lipidomics 

##Objective: 
#1. Check to make sure all lipid has the correct number of side chains identified
#2. Remove duplicated lipids

##Steps:
#1. Make sure all lipids have the right number of side chains, result: all lipid have the correct number of sidechains 
#2. Get mean intensity of each lipid across all samples 
#3. Remove duplicated lipids (lipids that have the identical headgroup and side chain identification) by only keeping one lipid that has the highest average intensity 


rm(list = ls())
library(tidyverse)
library(stringi)
library(ggpubr)
setwd(rstudioapi::getActiveProject())

load("./Output_Data/M2PM_pos_ion_clean.Rdata") #1122
load("./Output_Data/M2PM_neg_ion_clean.Rdata") #496

Ion.Clean.all <- bind_rows(M2PM.pos.clean, M2PM.neg.clean)

Cla.ls <- unique(Ion.Clean.all$Class)
one.chain.cla <- c("LPE", "LPI", "LPC", "Cholesterol", "AcCa", "MG", "Co", "LPS")
two.chain.cla <- c("PE", "PI", "PC", "Cer", "SM", "PS", "PG", "DG")
tri.chain.cla <- "TG"
quad.chain.cla <- "CL"

## For each lipid class, only keep lipids that have expected number of side chains detected
chain.cln <- Ion.Clean.all %>% 
  rowwise() %>% 
  mutate(Chain.sep = str_count(LipidIon, "\\_")) %>% 
  mutate(chain.check = case_when(
    Class %in% one.chain.cla & Chain.sep == 0 ~ "T",
    Class %in% two.chain.cla & Chain.sep == 1 ~ "T",
    Class %in% tri.chain.cla & Chain.sep == 2 ~ "T",
    Class %in% quad.chain.cla & Chain.sep == 3 ~ "T",
  )) %>% 
  relocate(c("Chain.sep","chain.check"), .after = LipidIon) %>%
  filter(chain.check == "T") %>% #426 lipids total
  mutate(LipidnoRT = LipidIon) 

#================#================#================#================
## Remove duplicated lipids, 426 ion cleaned lipids, only 245 of them have unique annotation
ls <- as.data.frame(table(chain.cln$LipidIon))

dup.cln <- chain.cln %>% #181 entry with duplicated lipidion ID
  filter(LipidIon %in% as.character(ls$Var1[ls$Freq>1]))  #181

dup.cln$AvgInt.smpl <- rowMeans(subset(dup.cln, select = XZ_43:XZ_82), na.rm = TRUE) #for each lipid with identical ID, pick the one with the highest overall intensity across all **experimental samples**

dup.cln <- dup.cln %>% #78, after removing 103 duplicated entries
  group_by(LipidIon) %>% 
  group_modify( ~{
    .x %>% 
      filter(AvgInt.smpl == max(AvgInt.smpl))
  })

uni.cln <- chain.cln %>% 
  filter(LipidIon %in% as.character(ls$Var1[ls$Freq==1])) #lipids with unique entry = 245

M2PM.cln.all <- bind_rows(uni.cln, dup.cln) #323 lipids after final removal of duplicated annotation

save(M2PM.cln.all, file = "./Output_Data/M2PM_LipMol_Ion+FA_clean_323.Rdata")
