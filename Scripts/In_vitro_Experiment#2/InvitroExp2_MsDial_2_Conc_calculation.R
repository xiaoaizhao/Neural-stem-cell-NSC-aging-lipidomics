### Concentration calculation

setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)

Exp2 <- read.csv("./Input_Data/Exp2_Reintegration_values_121725.peak.color.coded.csv", stringsAsFactors = F)
load("./Output_Data/Exp2.MS-Dial_IS.all.Rdata")

Exp2.rmv.poor.peaks <- Exp2 %>% 
  filter(!Peak.color.code == "Orange") #278

E2.lpd.check <- Exp2.rmv.poor.peaks %>% 
  rowwise() %>% 
  mutate(Class = str_split(Metabolite.name, " ")[[1]][1]) %>% 
  group_by(Metabolite.name) %>% 
  group_modify(~{
    .x %>% 
      mutate(Ion.check = ifelse(Adduct.type == IS2.export.flt$Adduct.type[IS2.export.flt$Class == Class], "T", "F"))  
  }) %>% 
  relocate(c(Class, Ion.check), .after = "Metabolite.name") %>% 
  filter(Ion.check == "T") ##all 278 lipids checked out

## transform into long form
Lpd.exp2 <- E2.lpd.check %>% 
  ungroup() %>% 
  select(c(Metabolite.name, Metabolite, Adduct.type, matches("^O|^Y"))) %>% 
  select(-Ontology) %>% 
  pivot_longer(-c(Metabolite.name, Adduct.type, Metabolite), names_to = "Sample", values_to = "Intensity") %>% 
  rowwise() %>% 
  mutate(Class = str_split(Metabolite.name, " ")[[1]][1])

##append IS value of each sample with corresponding spike in concentration
##single point concentration calculation = [detected_lipid_intensity] /([internal_standard_intensity]/[concentration_of_internal_standard])
load("./Output_data/IS_Exp2_w_conc_MsD.Rdata")
df.conc <- Lpd.exp2 %>%
  mutate(., smpl = Sample) %>%
  mutate(., Class1 = Class) %>%
  mutate(., Ion1 = Adduct.type) %>%
  group_by(Class, Adduct.type, Sample) %>%
  group_modify(~ {
    .x %>%
      mutate(., IS_int = as.numeric(IS_Exp2_w_conc_MsD[IS_Exp2_w_conc_MsD$Class == unique(Class1) 
                                                   & IS_Exp2_w_conc_MsD$Adduct.type == unique(Ion1), 
                                                   unique(smpl)])) %>%
      mutate(., conc_in_smpl = as.numeric(IS_Exp2_w_conc_MsD[IS_Exp2_w_conc_MsD$Class == unique(Class1)
                                                         & IS_Exp2_w_conc_MsD$Adduct.type  == unique(Ion1), "Conc_in_sample"])) %>%
      mutate(., Endo_conc = Intensity/(IS_int/conc_in_smpl))
  })


### sanity check
s.n <- "Y4_M"
Cla.n <- "PE"
r.ind <- 10
a<- df.conc$Intensity[df.conc$Class == Cla.n & df.conc$Sample == s.n][r.ind]/df.conc$Endo_conc[df.conc$Class == Cla.n & df.conc$Sample == s.n][r.ind]
b <- as.numeric(IS_Exp2_w_conc_MsD[IS_Exp2_w_conc_MsD$Class == Cla.n, s.n])/IS_Exp2_w_conc_MsD$Conc_in_sample[IS_Exp2_w_conc_MsD$Class == Cla.n]
round(a,6) == round(b,6)

df.conc.wide <- df.conc %>%
  ungroup() %>% 
  select(., c("Sample", "Metabolite.name", "Metabolite", "Endo_conc")) %>%
  pivot_wider(names_from = Sample, values_from = Endo_conc)

Exp2.lipid.w.conc.MsD <- df.conc.wide
save(Exp2.lipid.w.conc.MsD, file = paste0("./Output_data/Exp2_278.quant.lipids.byMsD.Rdata"))

