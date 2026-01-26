### Extract internal standards from all samples of Mboat2 overexpression and plasma membrane lipid supplementation lipidomics

setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)
library(matrixStats)

##Input IS concentration 

## ==== Ion list for internal standards ====
load("./Output_Data/IS_Ion_for_quantification.Rdata")

## ==== Mboat2 OE and lipid supplementation ====
## ==== Mboat2 OE and lipid supplementation ====

M2OE.label <- read.csv("./Input_Data/M2OE_validated_standards[63].csv", stringsAsFactors = F)

Smp.Key <- read.csv("./Input_Data/March_Sample_list_071123_forR.csv", stringsAsFactors = F)

IS.M2OE<- M2OE.label %>% 
  mutate(IS.ion = paste0(Metabolite.name, Adduct.type)) %>% 
  filter(IS.ion %in% Ion.list.per.class ) %>% 
  select(Metabolite.name, Adduct.type, IS.ion, starts_with("XZ_")) %>% 
  rename_at(vars(matches(Smp.Key$Sample_ID)), ~Smp.Key$Sample.Name)

save(IS.M2OE, file = "./Output_Data/M2OE.MS-Dial_IS.all.Rdata")
