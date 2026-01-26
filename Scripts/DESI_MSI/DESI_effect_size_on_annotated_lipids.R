### Match annotated DESI lipids with In vitro and In vivo LC-MS/MS datasets

rm(list = ls())
library(tidyverse)
library(stringi)
source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")
##====match In vitro ====


Invitro.new.peak <- read.csv("./Input_Data/Exp3_Validated_Reintegration_121125_color.code.csv", stringsAsFactors = F)

Invitro <- Invitro.new.peak %>% 
  filter(!Peak.color.code == "orange") %>% #404 lipids
  rowwise() %>% 
  mutate(ion = substr(Adduct.type, str_locate(Adduct.type, "M")+1, str_locate(Adduct.type, "]")-1)) %>% 
  rename("LipidIon" = "Metabolite.name")

Invitro.n <- MsD.lpd.rmv.abc(Invitro)

##====match in vivo====
Invivo.df <- read.csv("./Input_Data/Invivo_sorted_NSCs_integration_graded_120325.csv", stringsAsFactors = F)
Invivo.peak <- Invivo.df %>% 
  filter(Peak_shape == "Yes") #29

Invivo.df <- Invivo.peak %>% 
  mutate(Metabolite.name = ifelse(Alignment.ID == 3476, "ST 27:1;O", Metabolite.name)) %>% 
  mutate(Metabolite = ifelse(Alignment.ID == 3476, "ST 27:1;O", Metabolite)) %>% 
  rowwise() %>% 
  mutate(ion = substr(Adduct.type, str_locate(Adduct.type, "M")+1, str_locate(Adduct.type, "]")-1)) %>% 
  rename("LipidIon" = "Metabolite.name")

Invivo.n <- MsD.lpd.rmv.abc(Invivo.df)
 
##====Take the average of m.z from 3 studies and convert to -H for overlapping with DESI====
mz.2data <- bind_rows(Invitro.n, Invivo.n) %>% 
  group_by(LipidIon, ion) %>%
  summarise(MZ_avg = mean(Average.Mz, na.rm = TRUE))

mz.2data.negH <- Conv.NegH(mz.2data)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DESI.ls <- read.csv("./Input_Data/DESI_validated_lipids_list_072922.csv", stringsAsFactors = F)

DESI.ls <- DESI.ls %>% 
  rename(LipidIon = Ion)

DESI.rename <- ether.rename(DESI.ls)
ID.list <- DESI.rename$LipidIon
ID.list[ID.list == "SM(d18:1_18:0)"] <- "SM(18:1;O2/18:0)"

lpd.anno.inDESI <- mz.2data.negH %>% 
  filter(., LipidIon %in% ID.list) # 11

## Read in effect size matrix of DESI dataset
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DESI <- read.csv("./Output_Data/20210512_DESI_decomposition_ES_oldVsYoung.csv", stringsAsFactors = F)

DESI_GFAP <- DESI %>%
  filter(., cell == "gfap") %>%
  dplyr::select(c("peak", "es", "se" ,"pv"))

ES.lpd.anno.inDESI <- inner_join(DESI_GFAP, lpd.anno.inDESI, by = "peak") 
save(ES.lpd.anno.inDESI, file = "./Output_Data/EffectSize_DESI_annotated_lipids.Rdata")

