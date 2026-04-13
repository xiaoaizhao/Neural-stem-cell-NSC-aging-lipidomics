### Extract internal standards from all samples of _In vitro_ lipidomics

## ------------------------------------------------------------------------------------------------------
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
Neg3 <- read.csv("./Input_Data/Neg_Cust_invitro.csv", stringsAsFactors = F)
Pos3 <- read.csv("./Input_Data/Pos_Cust_invitro.csv", stringsAsFactors = F)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
PI <- read.csv("./Input_Data//Exp3_PI_IS_ID6845.csv", stringsAsFactors = F) %>% 
  filter(Alignment.ID == 6845)
Cholesterol <- read.csv("./Input_Data/Exp3_Cholesterol_IS_ID2956.csv", stringsAsFactors = F) %>% 
  filter(Alignment.ID == 2956)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
Neginvitro <- bind_rows(Neg3 %>% filter(!grepl("PI", Metabolite.name)), PI) 
Posinvitro <- bind_rows(Pos3, Cholesterol) %>% 
  mutate(Adduct.type = ifelse(grepl("Cer", Metabolite.name), "[M+H-H2O]+", Adduct.type)) 


## ---------------------------------------------------------------------------------------------------------------------------------------------------------

IS.invitro <- bind_rows(Posinvitro, Neginvitro) %>% 
  filter(grepl("Cer|LPE|PE|SM", Metabolite.name)) %>% 
  select(c(Metabolite.name, Adduct.type), matches("^QC|^XZ")) %>% 
  select(-contains("MSMS"))
IS.invitro.L <- IS.invitro %>%
  pivot_longer(-c(Metabolite.name, Adduct.type), names_to = "Sample", values_to = "Int") %>% 
  mutate(Mode = ifelse(substr(Adduct.type, nchar(Adduct.type), nchar(Adduct.type)) == "+", "Pos", "Neg"))

IS.invitro.pos <- IS.invitro.L %>% 
  filter(Mode == "Pos") %>% 
  rename("Int.pos" = "Int")

IS.invitro.neg <- IS.invitro.L %>% 
  filter(Mode == "Neg") %>% 
  rename("Int.neg" = "Int")

IS.invitro.all <- inner_join(IS.invitro.pos, IS.invitro.neg, by = c("Metabolite.name", "Sample")) %>% 
  mutate(Log2intpos = log2(Int.pos)) %>% 
  mutate(Log2intneg = log2(Int.neg)) 


## ---------------------------------------------------------------------------------------------------------------------------------------------------------

mean.int.3 <- IS.invitro.L %>% 
  group_by(Metabolite.name, Mode) %>% 
  summarise(MeanInt = mean(Int))

IS3.adduct <- IS.invitro.L %>%
  select(Metabolite.name, Adduct.type, Mode) %>% 
  group_by(Metabolite.name, Mode) %>% 
  summarise(Adduct = unique(Adduct.type))

mean.int.3<- left_join(mean.int.3, IS3.adduct, by = c("Metabolite.name", "Mode"))
mean.int.3max.int <- mean.int.3 %>% 
  group_by(Metabolite.name) %>% 
  filter(! (Metabolite.name == "Cer_NS 18:1(d7);2O/15:0" & Adduct == "[M+HCOO]-")) %>% 
  filter(! (Metabolite.name == "LPE 18:1(d7)" & Adduct == "[M+H]+")) %>% 
  filter(! (Metabolite.name == "PE 15:0_18:1(d7)" & Adduct == "[M-H]-")) %>% 
  filter(! (Metabolite.name == "SM 18:1;2O/18:1(d9)" & Adduct == "[M+HCOO]-"))

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
IS.invitro.export <-bind_rows(Posinvitro, Neginvitro)

IS.invitro.export.flt <- left_join(IS.invitro.export, mean.int.3max.int, by = c("Metabolite.name")) %>% 
  mutate(Filter = case_when(
    is.na(Adduct) ~ "Keep",
    !is.na(Adduct)&Adduct == Adduct.type ~ "Keep"
  )) %>% 
  filter(Filter == "Keep") %>% 
  select(c(Metabolite.name, Adduct.type, Adduct), matches("^QC|^XZ")) %>% 
  mutate(Class = substr(Metabolite.name, 1, str_locate(Metabolite.name, " ")-1)) %>% 
  mutate(Class = ifelse(grepl("^Cer", Metabolite.name), "Cer", Class)) %>% 
  mutate(Class = ifelse(grepl("^ST", Metabolite.name), "Cholesterol", Class)) %>% 
  relocate(Class, .after = "Metabolite.name")  %>% 
  select(-contains("MSMS")) %>% 
  filter(!grepl("CL", Metabolite.name))

save(IS.invitro.export.flt, file = "./Output_Data/Invitro.MS-Dial_IS.all.Rdata")

## Export ions for each IS for quantitative analysis
load("./Output_Data/Invitro.MS-Dial_IS.all.Rdata")
IS.ion.quant <- IS.invitro.export.flt %>% 
  mutate(IS.ion = paste0(Metabolite.name, Adduct.type))

Ion.list.per.class <- IS.ion.quant$IS.ion
save(Ion.list.per.class, file = "./Output_Data/IS_Ion_for_quantification.Rdata")
