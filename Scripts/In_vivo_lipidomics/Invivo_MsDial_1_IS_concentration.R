### Append concentration to each internal standards

rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

load("./Output_Data/IS_invivo_samples_Int.MsD.Rdata")

IS.label <- read.csv("./Input_Data/Labelled_IS_concentration_2020_052020.csv", stringsAsFactors = F,
                     check.names = F, colClasses = c("character", rep("numeric", 4)))

Label.conc <- IS.label %>%
  mutate(., Class = c("PC", "LPC", "PE", "LPE", "PG", "PI", "PS", "TG", "DG", "MG", "ChE", "SM", "Cer", "Cholesterol")) %>%
  mutate(., dilution = `Conc. (uM)` / `Conc_in_sample`)
##add labelled lipid concentration to IS matrix
Invivo.IS.Conc <- left_join(Invivo.sample.IS, Label.conc, by="Class") %>% 
  select(-c("Mixture Component","Molecular Weight","Conc. (ug/mL)","dilution")) %>% 
  relocate(Conc_in_sample, .before = "Conc. (uM)")

save(Invivo.IS.Conc, file = "./Output_Data/Invivo_IS_w_conc.Rdata")

