### Extract internal standards from all samples of _In vitro_ Experiment #2 lipidomics
##------------------------------------------------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(ggthemes)
library("scales")
library(ggpubr)
library(ggbeeswarm)
library(stringi)
library(openxlsx)
setwd(rstudioapi::getActiveProject())

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
Neg2 <- read.csv("./Input_Data/Neg_Cust_2.csv", stringsAsFactors = F)
Pos2 <- read.csv("./Input_Data/Pos_Cust_2.csv", stringsAsFactors = F)
DG <- read.csv("./Input_Data/DG-d7(manual peak integration).csv", stringsAsFactors = F)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
Pos2.fix.p <- Pos2 %>% 
  mutate(Adduct.type = ifelse(grepl("Cer", Metabolite.name), "[M+H-H2O]+", Adduct.type)) %>% 
  filter(!grepl("DG", Metabolite.name))
IS2 <- bind_rows(Pos2.fix.p,DG, Neg2) %>% 
  filter(grepl("Cer|LPE|PE|SM", Metabolite.name)) %>% 
  select(c(Metabolite.name, Adduct.type), matches("^Y|^O")) %>% 
  select(-Ontology)
IS2L <- IS2 %>%
  pivot_longer(-c(Metabolite.name, Adduct.type), names_to = "Sample", values_to = "Int") %>% 
  mutate(Mode = ifelse(substr(Adduct.type, nchar(Adduct.type), nchar(Adduct.type)) == "+", "Pos", "Neg"))

IS2.pos <- IS2L %>% 
  filter(Mode == "Pos") %>% 
  rename("Int.pos" = "Int")

IS2.neg <- IS2L %>% 
  filter(Mode == "Neg") %>% 
  rename("Int.neg" = "Int")

IS2.all <- inner_join(IS2.pos, IS2.neg, by = c("Metabolite.name", "Sample")) %>% 
  mutate(Log2intpos = log2(Int.pos)) %>% 
  mutate(Log2intneg = log2(Int.neg)) 


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
mean.int.2 <- IS2L %>% 
  group_by(Metabolite.name, Mode) %>% 
  summarise(MeanInt = mean(Int))

IS2.adduct <- IS2L %>%
  select(Metabolite.name, Adduct.type, Mode) %>% 
  group_by(Metabolite.name, Mode) %>% 
  summarise(Adduct = unique(Adduct.type))

mean.int.2<- left_join(mean.int.2, IS2.adduct, by = c("Metabolite.name", "Mode"))
mean.int.2max.int <- mean.int.2 %>% 
  group_by(Metabolite.name) %>% 
  group_modify( ~{
    .x %>% 
      filter(MeanInt == max(MeanInt))
  })


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
IS2.export <- bind_rows(Pos2.fix.p,DG, Neg2)

IS2.export.flt <- left_join(IS2.export, mean.int.2max.int, by = c("Metabolite.name")) %>% 
  mutate(Filter = case_when(
    is.na(Adduct) ~ "Keep",
    !is.na(Adduct)&Adduct == Adduct.type ~ "Keep"
  )) %>% 
  filter(Filter == "Keep") %>% 
  select(c(Metabolite.name, Adduct.type, Adduct), matches("^Y|^O|^QC")) %>% 
  mutate(Class = substr(Metabolite.name, 1, str_locate(Metabolite.name, " ")-1)) %>% 
  mutate(Class = ifelse(grepl("^Cer", Metabolite.name), "Cer", Class)) %>% 
  relocate(Class, .after = "Metabolite.name") %>% 
  select(-Ontology) %>% 
  filter(!grepl("CL", Metabolite.name))

save(IS2.export.flt, file = "~/Dropbox/NSC_rebuttal/Output_Data/Exp2.MS-Dial_IS.all.Rdata")

