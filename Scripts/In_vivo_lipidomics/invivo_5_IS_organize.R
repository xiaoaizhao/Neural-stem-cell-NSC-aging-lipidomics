
##Extract Internal Standard for each lipid class
##This script append IS concentration in sample to the matrix
##0.5μl of EquiSplash + 0.1μl cholesterol were used in each sample, final reconsitituted volume = 55μl
##For lipid species included in EquiSplash there is a 110 fold dilution, Cholesterol standard has 550 fold dilution, identical to 2020 GPMV run 
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

##match IS reading with sample name
load(file = "./Output_Data/IS_all_sample.Rdata")

label <- read.csv("./Input_Data/Lable_sample_ID.csv", stringsAsFactors = F)

All.IS.rename <- All.IS %>% 
  rename_at(vars(label$Running_ID[1:12]), ~label$Sorted.cell.samples[1:12])


IS.df <- All.IS.rename %>%
  mutate(., Class = c("Cer", "Cholesterol", "ChE", "DG", "LPC", "LPE", "PC", "PE", "PG", "PI", "PS", "SM", "TG")) %>%
  mutate(., Ion = substr(IS, nchar(IS)-1, nchar(IS))) %>%
  mutate(., Ion = ifelse(grepl("+NH4", IS), "+NH4", Ion)) %>%
  mutate(., Ion = ifelse(grepl("+H-H2O", IS), "+H-H2O", Ion))

##Add "[]" in ion to match data matrix
IS.df <- IS.df %>%
  mutate(., Ion = paste0("[", Ion, "]"))

##load in labelled lipid concentration
IS.label <- read.csv("./Input_Data/Labelled_IS_concentration_2020_052020.csv", stringsAsFactors = F,
                     check.names = F, colClasses = c("character", rep("numeric", 4)))

Label.conc <- IS.label %>%
  mutate(., Class = c("PC", "LPC", "PE", "LPE", "PG", "PI", "PS", "TG", "DG", "MG", "ChE", "SM", "Cer", "Cholesterol")) %>%
  mutate(., dilution = `Conc. (uM)` / `Conc_in_sample`)
##add labelled lipid concentration to IS matrix
IS_w_conc <- left_join(IS.df, Label.conc, by="Class")

save(IS_w_conc, file = paste0("./Output_Data/Invivo_13_IS_int_w_conc.Rdata"))

