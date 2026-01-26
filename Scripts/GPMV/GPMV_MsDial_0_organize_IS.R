### Extract internal standards from all samples of GPMV lipidomics

setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)
library(matrixStats)

##Input IS concentration 
load("./Output_Data/IS_Ion_for_quantification.Rdata")

## ==== GPMV ====
## ==== GPMV ====

GPMV.label <- read.csv("./Input_Data/Lipid_Standards_GPMV_112125.csv", stringsAsFactors = F)
GPMV.ls <- read.csv("./Input_Data/GPMV_sample_list.csv", stringsAsFactors = F)

GPMV.ls <- GPMV.ls %>% 
  filter(!grepl("XZ_3|XZ_10|XZ_18", ID))

IS.GPMV<- GPMV.label %>% 
  mutate(IS.ion = paste0(Metabolite.name, Adduct.type)) %>% 
  filter(IS.ion %in% Ion.list.per.class ) %>% 
  select(Metabolite.name, Adduct.type, IS.ion, starts_with("X")) %>% 
  select(-c("X1.1", "X1.2")) %>% 
  pivot_longer(-c(Metabolite.name, Adduct.type, IS.ion), names_to = "Sample", values_to = "Int") %>% 
  rowwise() %>% 
  mutate(Sample_name = paste0("XZ_", substr(Sample, 2, nchar(Sample)))) %>% 
  select(-Sample) %>% 
  pivot_wider(c(Metabolite.name, Adduct.type, IS.ion), names_from = "Sample_name", values_from = "Int") %>% 
  rename_at(vars(matches(GPMV.ls$ID)), ~GPMV.ls$Sample.list)

save(IS.GPMV, file = "./Output_Data/GPMV.MS-Dial_IS.all.Rdata")
