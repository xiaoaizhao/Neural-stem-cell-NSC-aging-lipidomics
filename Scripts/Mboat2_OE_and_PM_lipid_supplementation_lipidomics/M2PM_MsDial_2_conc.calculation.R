### Concentration calculation

#Mboat2 overexpression and lipid supplementation data
setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)

M2PMdf <- read.csv("./Input_Data/M2OE_Validation_112225.csv", stringsAsFactors = F)
load("./Output_Data/M2PM.IS.MsD.w.conc.Rdata")

M2PM.lpd.check <- M2PMdf %>% 
  rowwise() %>% 
  mutate(Class = str_split(Metabolite.name, " ")[[1]][1]) %>%
  mutate(Class = ifelse(Class == "ST", "Cholesterol", Class)) %>% 
  mutate(Class = ifelse(grepl("Cer", Class),"Cer", Class)) %>% 
  relocate(Class, .after = Metabolite.name) %>% 
  group_by(Metabolite.name) %>% 
  group_modify(~{
    .x %>% 
      mutate(Ion.check = ifelse(Adduct.type == M2PM.IS.MsD.conc$Adduct.type[M2PM.IS.MsD.conc$Class == Class], "T", "F"))  
  }) %>% 
  relocate(c(Class, Ion.check), .after = "Metabolite.name") %>% 
  filter(Ion.check == "T") ##all 397 lipids check out


## Append sample name to file name
Smp.Key <- read.csv("./Input_Data/March_Sample_list_071123_forR.csv", stringsAsFactors = F)

M2PM.lpd.int<- M2PM.lpd.check %>% 
  select(Metabolite.name, Adduct.type, starts_with("XZ_"), Class) %>% 
  pivot_longer(-c(Metabolite.name, Adduct.type, Class), names_to = "Sample", values_to = "Intensity") %>% 
  rowwise() %>% 
  mutate(Sample_ID = Smp.Key$Sample.Name[Smp.Key$Sample_ID == Sample]) %>% 
  select(-c(Sample)) %>% 
  rename("Sample" = "Sample_ID")

##append IS value of each sample with corresponding spike in concentration
##single point concentration calculation = [detected_lipid_intensity] /([internal_standard_intensity]/[concentration_of_internal_standard])

df.conc <- M2PM.lpd.int %>%
  mutate(., smpl = Sample) %>%
  mutate(., Class1 = Class) %>%
  mutate(., Ion1 = Adduct.type) %>%
  group_by(Class, Adduct.type, Sample) %>%
  group_modify(~ {
    .x %>%
      mutate(., IS_int = as.numeric(M2PM.IS.MsD.conc[M2PM.IS.MsD.conc$Class == unique(Class1) 
                                                   & M2PM.IS.MsD.conc$Adduct.type == unique(Ion1), 
                                                   unique(smpl)])) %>%
      mutate(., conc_in_smpl = as.numeric(M2PM.IS.MsD.conc[M2PM.IS.MsD.conc$Class == unique(Class1) 
                                                         & M2PM.IS.MsD.conc$Adduct.type  == unique(Ion1), "Conc_in_23_2Batches"])) %>%
      mutate(., Endo_conc = Intensity/(IS_int/conc_in_smpl))
  })


### sanity check
s.n <- "Y5_FLG_Mb2_OE"
Cla.n <- "LPC"
r.ind <- 3
a<- df.conc$Intensity[df.conc$Class == Cla.n & df.conc$Sample == s.n][r.ind]/df.conc$Endo_conc[df.conc$Class == Cla.n & df.conc$Sample == s.n][r.ind]
b <- as.numeric(M2PM.IS.MsD.conc[M2PM.IS.MsD.conc$Class == Cla.n, s.n])/M2PM.IS.MsD.conc$Conc_in_23_2Batches[M2PM.IS.MsD.conc$Class == Cla.n]
a == b

df.conc.wide <- df.conc %>%
  ungroup() %>% 
  select(., c("Sample", "Metabolite.name", "Endo_conc")) %>%
  pivot_wider(names_from = Sample, values_from = Endo_conc)

M2PM.lipid.w.conc.MsD <- df.conc.wide
save(M2PM.lipid.w.conc.MsD, file = paste0("./Output_data/M2PM_397.quant.lipids.byMsD.Rdata"))
