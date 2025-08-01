## Match annotated DESI lipids (by a separate tandem MS on the same set of samples) with their respective effect size by accurate mass (m/z)

## Steps:
#1. Calculate the average m/z in lipids from all 3 datasets
#2. Convert all identified lipids from Primary culture #1 - #3 to [-H] ion addut (using the average m/z from step #1)
#3. Match m/z from in vitro summary to m/z detected from DESI-MSI in [-H] ion adduct

rm(list = ls())
library(tidyverse)
library(stringi)
source("./Scripts/Function_scripts/Effect_size_functions.R")

##====match primary culture #1====
load(file = "./Output_data/2_Ion+FA_clean_373_lipids.Rdata")
exp1.df <- FA_unique373  %>% 
  ungroup() %>% 
  select(LipidIon, m.z, Ion) %>% 
  mutate(ion= substr(Ion, 2, str_locate(Ion, "\\]")-1)) %>% 
  mutate(Exp = "Exp1") %>% 
  rowwise() %>% 
  mutate(Ion.check = ifelse(grepl(ion, LipidIon), "T", "F")) %>% 
  filter(Ion.check == "T")

length(rownames(exp1.df)) == length(rownames(FA_unique373))

exp1.n <- ether.rename(exp1.df)


##====match primary culture #2====
load("./Output_Data/Exp2_Dup_ID_rmv+dup_isomer_rmv_694_lipids.Rdata")
exp2.df <- df.rmv.isomer %>% 
  ungroup() %>% 
  select(LipidIon, CalcMz, Ion) %>% 
  mutate(ion= substr(Ion, 2, str_locate(Ion, "\\]")-1)) %>% 
  mutate(Exp = "Exp2") %>% 
  rowwise() %>% 
  mutate(Ion.check = ifelse(grepl(ion, LipidIon), "T", "F")) %>% 
  filter(Ion.check == "T") %>% 
  rename("m.z" = "CalcMz")

length(rownames(exp2.df)) == length(rownames(df.rmv.isomer))

exp2.n <- ether.rename(exp2.df)

##====match primary culture #3====
load("./Output_Data/Exp3_LipMol_Ion+FA_clean_414.Rdata") ##lipid data frame with m/z for each lipid species

exp3.df <- Exp3.cln.all %>% #414, need to manually fix cholesterol annotation
  ungroup() %>% 
  select(LipidIon, m.z, Ion) %>% 
  mutate(LipidIon = ifelse(grepl("Ch\\+H\\-H2O", LipidIon), "Cholesterol(0:0)+H-H2O", LipidIon)) %>% 
  mutate(ion= substr(Ion, 2, str_locate(Ion, "\\]")-1)) %>% 
  mutate(Exp = "Exp3") %>% 
  rowwise() %>% 
  mutate(Ion.check = ifelse(grepl(ion, LipidIon), "T", "F")) %>% 
  filter(Ion.check == "T")

length(rownames(exp3.df)) == length(rownames(Exp3.cln.all))

exp3.n <- ether.rename(exp3.df)
  
##====Take the average of m.z from 3 studies and convert to -H for overlapping with DESI====
mz.3data <- bind_rows(exp1.n, exp2.n, exp3.n) %>% 
  group_by(LipidIon, ion) %>%
  summarise(MZ_avg = mean(m.z, na.rm = TRUE))

mz.3data.negH <- Conv.NegH(mz.3data)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DESI.ls <- read.csv("./Input_Data/DESI_validated_lipids_list_072922.csv", stringsAsFactors = F)

DESI.ls <- DESI.ls %>% 
  rename(LipidIon = Ion)

DESI.rename <- ether.rename(DESI.ls)
ID.list <- DESI.rename$LipidIon


lpd.anno.inDESI <- mz.3data.negH %>% 
  filter(., LipidIon %in% ID.list) # 13 out of 14 lipids from structurally verified DESI data were found in the effect size matrix

## Read in effect size matrix of DESI dataset
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DESI <- read.csv("./Output_Data/20210512_DESI_decomposition_ES_oldVsYoung.csv", stringsAsFactors = F)

DESI_GFAP <- DESI %>%
  filter(., cell == "gfap") %>%
  dplyr::select(c("peak", "es", "se" ,"pv"))

ES.lpd.anno.inDESI <- inner_join(DESI_GFAP, lpd.anno.inDESI, by = "peak") 
save(ES.lpd.anno.inDESI, file = "./Output_Data/EffectSize_DESI_annotated_lipids.Rdata")
