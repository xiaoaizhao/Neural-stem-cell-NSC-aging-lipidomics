### Concentration calculation

# GPMV data
setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)

GPMVdf <- read.csv("./Input_Data/GPMV_validated_120525_Full_integration_peak_graded.csv", stringsAsFactors = F)
load("./Output_Data/GPMV.IS.w.conc.Rdata")

GPMV.lpd.check <- GPMVdf %>% 
  rowwise() %>% 
  mutate(Class = str_split(Metabolite.name, " ")[[1]][1]) %>%
  mutate(Class = ifelse(Class == "ST", "Cholesterol", Class)) %>% 
  mutate(Class = ifelse(Class == "HexCer", "Cer", Class)) %>% 
  relocate(Class, .after = Metabolite.name) %>% 
  group_by(Metabolite.name) %>% 
  group_modify(~{
    .x %>% 
      mutate(Ion.check = ifelse(Adduct.type == GPMV.IS.w.conc$Adduct.type[GPMV.IS.w.conc$Class == Class], "T", "F"))  
  }) %>% 
  relocate(c(Class, Ion.check), .after = "Metabolite.name") %>% 
  filter(Ion.check == "T") ##all 332 lipids check out

# remove peaks with poor quality

GPMV.poor.peak.rmv <- GPMV.lpd.check %>% 
  filter(!Poor_peak_shape == "Yes") #291, removed 41 peaks with poor shapes

poor.peaks <- GPMV.lpd.check %>% 
  filter(Poor_peak_shape == "Yes")

## Append sample name to file name
GPMV.ls <- read.csv("./Input_Data/GPMV_sample_list.csv", stringsAsFactors = F)

GPMV.ls <- GPMV.ls %>% 
  filter(!grepl("XZ_3|XZ_10|XZ_18", ID))

GPMV.lpd.int<- GPMV.poor.peak.rmv %>% 
  select(Metabolite.name, Adduct.type, starts_with("X"), Class) %>% 
  select(-c("X1.1", "X1.2")) %>% 
  pivot_longer(-c(Metabolite.name, Adduct.type, Class), names_to = "Sample", values_to = "Intensity") %>% 
  rowwise() %>% 
  mutate(Sample_name = paste0("XZ_", substr(Sample, 2, nchar(Sample)))) %>% 
  mutate(Sample_ID = GPMV.ls$Sample.list[GPMV.ls$ID == Sample_name]) %>% 
  select(-c(Sample, Sample_name)) %>% 
  rename("Sample" = "Sample_ID")

##append IS value of each sample with corresponding spike in concentration
##single point concentration calculation = [detected_lipid_intensity] /([internal_standard_intensity]/[concentration_of_internal_standard])

df.conc <- GPMV.lpd.int %>%
  mutate(., smpl = Sample) %>%
  mutate(., Class1 = Class) %>%
  mutate(., Ion1 = Adduct.type) %>%
  group_by(Class, Adduct.type, Sample) %>%
  group_modify(~ {
    .x %>%
      mutate(., IS_int = as.numeric(GPMV.IS.w.conc[GPMV.IS.w.conc$Class == unique(Class1) 
                                                   & GPMV.IS.w.conc$Adduct.type == unique(Ion1), 
                                                   unique(smpl)])) %>%
      mutate(., conc_in_smpl = as.numeric(GPMV.IS.w.conc[GPMV.IS.w.conc$Class == unique(Class1) 
                                                         & GPMV.IS.w.conc$Adduct.type  == unique(Ion1), "Conc_in_sample"])) %>%
      mutate(., Endo_conc = Intensity/(IS_int/conc_in_smpl))
  })


### sanity check
s.n <- "Y1_1219"
Cla.n <- "PE"
r.ind <- 5
a<- df.conc$Intensity[df.conc$Class == Cla.n & df.conc$Sample == s.n][r.ind]/df.conc$Endo_conc[df.conc$Class == Cla.n & df.conc$Sample == s.n][r.ind]
b <- as.numeric(GPMV.IS.w.conc[GPMV.IS.w.conc$Class == Cla.n, s.n])/GPMV.IS.w.conc$Conc_in_sample[GPMV.IS.w.conc$Class == Cla.n]
a == b

df.conc.wide <- df.conc %>%
  ungroup() %>% 
  select(., c("Sample", "Metabolite.name", "Endo_conc")) %>%
  pivot_wider(names_from = Sample, values_from = Endo_conc)

GPMV.lipid.w.conc <- df.conc.wide
save(GPMV.lipid.w.conc, file = paste0("./Output_data/GPMV_291.quant.lipids.Rdata"))
