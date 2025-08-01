##Clean duplicated lipid from Primary culture #3 2023 lipidomics data

##Objective: 
#1. Check to make sure all lipid has the correct number of side chains identified, 524 lipids total
#2. Remove duplicated lipids

##Steps:
#1. Make sure all lipids have the right number of side chains, result: all lipid have the correct number of sidechains - 524 total
#2. Get mean intensity of each lipid across all samples 
#3. Remove duplicated lipids (lipids that have the identical headgroup and side chain identification) by only keeping one lipid that has the highest average intensity. Removed 110 lipids


rm(list = ls())
library(tidyverse)
library(stringi)
library(ggpubr)
setwd(rstudioapi::getActiveProject())

load("./Output_Data/Exp3_pos_ion_clean.Rdata")#1095
load("./Output_Data/Exp3_neg_ion_clean.Rdata")#305

Ion.Clean.all <- bind_rows(Exp3.pos.clean, Exp3.neg.clean) #1400

Cla.ls <- unique(Ion.Clean.all$Class)
one.chain.cla <- c("LPE", "LPI", "LPC", "ChE", "Cholesterol", "AcCa", "Co", "MG", "LPS")
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
  filter(chain.check == "T") %>% #524 lipids total
  mutate(LipidnoRT = LipidIon) 

#================#================#================#================
## Remove duplicated lipids, 524 ion cleaned lipids, only 315 of them have unique annotation
ls <- as.data.frame(table(chain.cln$LipidIon))

dup.cln <- chain.cln %>% #209 entry with duplicated lipidion ID
  filter(LipidIon %in% as.character(ls$Var1[ls$Freq>1]))  #209

dup.cln$AvgInt.smpl <- rowMeans(subset(dup.cln, select = XZ_83:XZ_104), na.rm = TRUE) #for each lipid with identical ID, pick the one with the highest overall intensity across all **experimental samples**

dup.cln <- dup.cln %>% #99, after removing 110 duplicated entries
  group_by(LipidIon) %>% 
  group_modify( ~{
    .x %>% 
      filter(AvgInt.smpl == max(AvgInt.smpl))
  })

uni.cln <- chain.cln %>% 
  filter(LipidIon %in% as.character(ls$Var1[ls$Freq==1])) #lipids with unique entry = 315

Exp3.cln.all <- bind_rows(uni.cln, dup.cln) #414 lipids after final removal of duplicated annotation

save(Exp3.cln.all, file = "./Output_Data/Exp3_LipMol_Ion+FA_clean_414.Rdata")
