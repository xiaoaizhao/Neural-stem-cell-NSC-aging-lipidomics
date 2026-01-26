
### Append concentration to each internal standards

## Mboat2 OE and lipid supplementation
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

load("./Output_Data/M2OE.MS-Dial_IS.all.Rdata")
load("./Input_Data/IS_new_conc_July_batch1_2.Rdata")
Class <- as_tibble(c("PC", "LPC", "PE", "LPE", "PG", "PI", "PS", "TG", "DG", "MG", "CE", "SM", "Cer", "Cholesterol"))

conc.M2PM <- bind_cols(IS.Jul23.Batch1_2, Class) %>% 
  rename("Class" = "value")

M2PM.MsD <- IS.M2OE %>% 
  rowwise() %>% 
  mutate(Class = str_split(Metabolite.name, " ")[[1]][1]) %>% 
  relocate(Class, .after = "Metabolite.name") %>% 
  mutate(Class = ifelse(grepl("^Cer", Class), "Cer", Class)) %>% 
  mutate(Class = ifelse(grepl("^ST", Class), "Cholesterol", Class)) %>% 
  select(-c(IS.ion, starts_with("QC")))

M2PM.IS.MsD.conc <- left_join(M2PM.MsD, conc.M2PM, by = "Class") %>% 
  select(-c("Metabolite.name","Mixture Component","Molecular Weight","Conc. (ug/mL)","Conc_in_sample")) %>% 
  relocate("Conc_in_23_2Batches", .before = "Conc. (uM)")

save(M2PM.IS.MsD.conc, file = "./Output_Data/M2PM.IS.MsD.w.conc.Rdata")


