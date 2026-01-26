### Extract internal standards from all samples of _In vivo_ lipidomics

setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)
library(matrixStats)

## ==== IS from in vivo sorted cells ====
load("./Output_Data/IS_Ion_for_quantification.Rdata")

## ==== For in vivo dataset use [M-H]- to quantify Ether PE. Note: no PE species were validated

Ion.list.per.class[4] <- str_replace_all(Ion.list.per.class[4], "\\+", "\\-")

IS.label <- read.csv("./Input_Data/Validated_standards_Invivo_sort_NSCs_Updated.csv", stringsAsFactors = F)

IS.invivo<- IS.label %>% 
  mutate(IS.ion = paste0(Metabolite.name, Adduct.type)) %>% 
  filter(IS.ion %in% Ion.list.per.class )

label <- read.csv("./Input_Data/Lable_sample_ID.csv", stringsAsFactors = F)
label.df <- label %>% 
  filter(!Running_ID == "XY.18")
IS.sample <- IS.invivo %>% 
  rename_with(~ str_remove(.x, "pos_")) %>% 
  select(c(Metabolite.name, matches("^XY.|QC"), Adduct.type)) %>% 
  rename_at(vars(matches(label.df$Running_ID)), ~label.df$Sorted.cell.samples) %>% 
  pivot_longer(-c(Metabolite.name, Adduct.type), names_to = "Sample", values_to = "IS_int") %>% 
  mutate(Cat = ifelse(grepl("^QC", Sample), "QC", "Experiment_sample"))

Invivo.only <- IS.sample %>% 
  filter(grepl("^O|^Y", Sample)) %>% 
  select(-Cat) %>% 
  rowwise() %>% 
  mutate(Class = str_split(Metabolite.name, " ")[[1]][1]) %>% 
  mutate(Class = ifelse(grepl("^Cer|^ST", Class), case_when(
    grepl("^Cer", Class) ~ "Cer",
    grepl("^ST", Class) ~ "Cholesterol"
  ), Class)) %>% 
  select(-Metabolite.name) %>% 
  rename("IS_MsD" = "IS_int")

Invivo.sample.IS <- Invivo.only %>% 
  pivot_wider(names_from = "Sample", values_from = IS_MsD)
save(Invivo.sample.IS, file = "./Output_Data/IS_invivo_samples_Int.MsD.Rdata")

