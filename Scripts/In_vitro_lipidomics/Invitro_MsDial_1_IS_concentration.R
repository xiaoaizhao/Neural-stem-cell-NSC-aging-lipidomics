### Append concentration to each internal standards
setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)

## Generate data frame with concentration for each internal standard

IS.label <- read.csv("./Input_Data/Labelled_IS_concentration_2020_052020.csv", stringsAsFactors = F,
                     check.names = F, colClasses = c("character", rep("numeric", 4)))

IS.invitro.conc<- IS.label %>% 
  mutate(., `Conc. (uM)` = ifelse(`Mixture Component` == "Cholesterol(d7)", 1304.53, `Conc. (uM)`)) %>% 
  mutate(., Conc_in_Batch3 = ifelse(
    `Mixture Component` == "Cholesterol(d7)", 
    `Conc. (uM)` *0.2/120,
    `Conc. (uM)` *1/120
  ))

save(IS.invitro.conc, file = "./Output_data/IS_new_conc_Invitro.Rdata")

### Now organize data matrix of all deuterated IS standards
rm(list = ls())

load("./Output_data/IS_new_conc_Invitro.Rdata")
load("./Output_Data/Invitro.MS-Dial_IS.all.Rdata") # from Script #0

IS.invitro.conc.c <- IS.invitro.conc %>% 
  mutate(Class = c("PC", "LPC", "PE", "LPE", "PG", "PI", "PS", "TG", "DG", "MG", "CE", "SM", "Cer", "Cholesterol")) %>% 
  relocate(Class, .after = `Mixture Component`)

IS_w_conc.Invitro <- left_join(IS.invitro.export.flt, IS.invitro.conc.c, by = "Class") %>% 
  relocate(`Mixture Component`, .after = Metabolite.name)

save(IS_w_conc.Invitro, file = paste0("./Output_data/IS_Invitro_w_conc_MsD.Rdata"))
