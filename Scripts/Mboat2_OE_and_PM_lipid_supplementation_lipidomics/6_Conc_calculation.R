### Processing script 6, _Mboat2_ overexpression and young plasma membrane lipid supplementation lipidomics 
#1. Calculate concentration (remove IS labelled from endogenous data) on 323 lipids (ion cleaned, duplicated annotation removed)

setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(reshape2)
library(stringi)

load("./Output_data/IS_M2PM_w_conc.Rdata")

load("./Output_data/M2PM_LipMol_Ion+FA_clean_323.Rdata")

M2PM.all <- M2PM.cln.all %>% 
  mutate(., Ion.m = substr(Ion, 2, nchar(Ion)-1))

lipid.no.conc <- M2PM.all %>%
  filter(., !Class %in% IS_w_conc.M2PM$Class) #48 lipids have no IS for quantification

lipid.conc <- M2PM.all %>%
  filter(., Class %in% IS_w_conc.M2PM$Class) #275 lipids have IS for quantification

lipidmelt <- lipid.conc %>%   ## 11825 rows
  ungroup() %>% 
  select(., c(Class, LipidIon, Ion.m, QC_1:QC_3, XZ_43:XZ_82)) %>% 
  pivot_longer(-c(LipidIon, Class, Ion.m), names_to = "Samples", values_to = "Intensity")


##append IS value of each sample with corresponding spike in concentration
##single point concentration calculation = [detected_lipid_intensity] /([internal_standard_intensity]/[concentration_of_internal_standard])

df.conc <- lipidmelt %>%
  mutate(., smpl = Samples) %>%
  mutate(., Class1 = Class) %>%
  mutate(., Ion1 = Ion.m) %>%
  group_by(Class, Ion.m, Samples) %>%
  group_modify(~ {
    .x %>%
      mutate(., IS_int = as.numeric(IS_w_conc.M2PM[IS_w_conc.M2PM$Class == unique(Class1) 
                                                    & IS_w_conc.M2PM$Ion == unique(Ion1), 
                                                    unique(smpl)])) %>%
      mutate(., conc_in_smpl = as.numeric(IS_w_conc.M2PM[IS_w_conc.M2PM$Class == unique(Class1) 
                                                          & IS_w_conc.M2PM$Ion == unique(Ion1), "Conc_in_23_2Batches"])) %>%
      mutate(., Endo_conc = Intensity/(IS_int/conc_in_smpl))
  })

### sanity check
s.n <- "XZ_49"
Cla.n <- "PE"
r.ind <- 18
a<- df.conc$Intensity[df.conc$Class == Cla.n & df.conc$Samples == s.n][r.ind]/df.conc$Endo_conc[df.conc$Class == Cla.n & df.conc$Samples == s.n][r.ind]
b <- as.numeric(IS_w_conc.M2PM[IS_w_conc.M2PM$Class == Cla.n, s.n])/IS_w_conc.M2PM$Conc_in_23_2Batches[IS_w_conc.M2PM$Class == Cla.n]
a == b

df.conc.wide <- df.conc %>%
  ungroup() %>% 
  select(., c("LipidIon", "Endo_conc", "Samples")) %>%
  pivot_wider(names_from = Samples, values_from = Endo_conc)

M2PM.lipid.w.conc <- df.conc.wide
save(M2PM.lipid.w.conc, file = paste0("./Output_data/M2PM_275_conc_lipids.Rdata")) 

lipid.no.conc.c <- lipid.no.conc %>% 
  ungroup() %>% 
  select(., c(LipidIon, QC_1:QC_3, XZ_43:XZ_82)) #48

M2PM_lpd.all <- bind_rows(M2PM.lipid.w.conc, lipid.no.conc.c)#323
save(M2PM_lpd.all, file = paste0("./Output_data/M2PM_275_conc_48_noconc_lipids.Rdata"))

