### Append concentration to each internal standards
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

##IS from previous analysis
load("./Output_Data/GPMV.MS-Dial_IS.all.Rdata")
IS.GPMV.all <- IS.GPMV %>% 
  rowwise() %>% 
  mutate(Class = str_split(Metabolite.name, " ")[[1]][1]) %>% 
  relocate(Class, .after = "Metabolite.name") %>% 
  mutate(Class = ifelse(grepl("^Cer", Class), "Cer", Class)) %>% 
  mutate(Class = ifelse(grepl("^ST", Class), "Cholesterol", Class))

##load in labelled lipid concentration
Label <- read.csv("./Input_Data/Labelled_IS_concentration_2020_052020.csv", stringsAsFactors = F,
                  check.names = F, colClasses = c("character", rep("numeric", 4)))

Label.conc <- Label %>%
  mutate(., Class = c("PC", "LPC", "PE", "LPE", "PG", "PI", "PS", "TG", "DG", "MG", "CE", "SM", "Cer", "Cholesterol"))
##add labelled lipid concentration to IS matrix
GPMV.IS.w.conc<- left_join(IS.GPMV.all, Label.conc, by="Class") %>% 
  select(-c("Metabolite.name", "IS.ion","Mixture Component","Molecular Weight","Conc. (ug/mL)")) %>% 
  relocate(Conc_in_sample, .before = "Conc. (uM)")

save(GPMV.IS.w.conc, file = "./Output_Data/GPMV.IS.w.conc.Rdata")

