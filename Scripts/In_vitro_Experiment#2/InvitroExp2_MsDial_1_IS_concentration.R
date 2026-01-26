### Append concentration to each internal standards

setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)

##Input IS concentration 

load("./Output_Data/Exp2.MS-Dial_IS.all.Rdata")

labelled.conc <- read.csv("./Input_Data/KO_Labelled_IS_concentration_082720.csv", stringsAsFactors = F,
                          check.names = F, colClasses = c("character", rep("numeric", 4)))

labelled.conc.tomerge <- labelled.conc %>%
  mutate(., Class = c("PC", "LPC", "PE", "LPE", "PG", "PI", 
                      "PS", "TG", "DG", "MG", "ChE", "SM", "Cer"))
IS_Exp2_w_conc_MsD <- left_join(IS2.export.flt, labelled.conc.tomerge, by= "Class") %>% 
  relocate(c(`Conc. (uM)`, Conc_in_sample), .after = "Class") %>% 
  select(-c("Mixture Component", "Molecular Weight",  "Conc. (ug/mL)"))
save(IS_Exp2_w_conc_MsD, file = paste0("./Output_data/IS_Exp2_w_conc_MsD.Rdata"))

