##Organized spike-in standard data frame for GPMV data

##Extract Internal Standard for each lipid class
##This script append IS concentration in sample to the matrix
##0.5μl of EquiSplash + 0.1μl cholesterol were used in each sample, final reconstituted volume = 55μl
##For lipid species included in EquiSplash there is a 110 fold dilution, Cholesterol standard has 550 fold dilution, identical to in vivo sorted NSC data 
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

##load in IS intensity detected in each sample
P_ISnew <- read.csv("./Input_Data/Progenesis_IS_updated_SM_052920.csv", stringsAsFactors = F)

P_ISnew <- P_ISnew[order(P_ISnew$LipidIon),]
GPMV.IS <- P_ISnew %>%
  mutate(., label = LipidIon) %>%
  mutate(., Class = c("DG", "PC", "PE","PG","PI", "PS","TG", "LPC", "LPE", "LPE", "Chol", "Cer", "Cer","SM") ) %>%
  mutate(., Ion = substr(label, nchar(label)-1, nchar(label))) %>%
  mutate(., Ion = ifelse(grepl("+NH4", label), "+NH4", Ion)) %>%
  mutate(., Ion = ifelse(grepl("+H-H2O", label), "+H-H2O", Ion))

##load in labelled lipid concentration
Label <- read.csv("./Input_Data/Labelled_IS_concentration_2020_052020.csv", stringsAsFactors = F,
                  check.names = F, colClasses = c("character", rep("numeric", 4)))

Label.conc <- Label %>%
  mutate(., Class = c("PC", "LPC", "PE", "LPE", "PG", "PI", "PS", "TG", "DG", "MG", "CE", "SM", "Cer", "Chol"))
##add labelled lipid concentration to IS matrix
GPMV.IS.w.conc<- left_join(GPMV.IS, Label.conc, by="Class")

save(GPMV.IS.w.conc, file = paste0("./Output_Data/GPMV_IS_org_w_conc.Rdata"))
