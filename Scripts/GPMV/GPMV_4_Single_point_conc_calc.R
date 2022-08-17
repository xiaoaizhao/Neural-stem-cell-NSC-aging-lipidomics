
##Quantify classes with standards to get concentration based on one spike-in per class

##Steps:
#1. Separate lipids belongs to classes with IS or not
#2. Calculate lipid concentration in classes of lipids with IS - 483 lipids (from 14 classes) 
#   Lipid concentration = [detected_lipid_intensity] /([internal_standard_intensity]/[concentration_of_internal_standard])
#3. Combine with the 82 lipids with intensity only, without concentration
rm(list=ls())
library(tidyverse)
library(reshape2)
setwd(rstudioapi::getActiveProject())

load(file = "./Output_Data/GPMV_Ion+FA_clean_565_lipids.Rdata") #565 lipids total
load(file = "./Output_Data/GPMV_IS_org_w_conc.Rdata") ## Internal Standard value from all sample, with corresponding spike in concentration

GPMV.IS.w.conc <- GPMV.IS.w.conc %>%
  mutate(., Ion = paste0("[", Ion, "]")) %>%
  mutate(., Class = ifelse(grepl("Chol", Class), "Cholesterol", Class))

lipid.df <- GPMV.FA_unique565 %>%
  select(., c("Ion","LipidIon", matches("Y._|O._"))) %>%
  select(., -contains("match"))

####Separate cleaned lipids into lipids with concentration and without concentration ####
lipid.no.conc <- lipid.df %>%
  filter(., !Class %in% GPMV.IS.w.conc$Class) #82 lipids have no IS for quantification

lipid.conc <- lipid.df %>%
  filter(., Class %in% GPMV.IS.w.conc$Class)  #483 lipids have corresponding IS for quantification

####Calculate lipid concentration in classes of lipids with IS####
lipidmelt <- melt(lipid.conc)
lipidmelt$variable <- as.character(lipidmelt$variable)

##single point concentration calculation = [detected_lipid_intensity] /([internal_standard_intensity]/[concentration_of_internal_standard])
df.conc <- lipidmelt %>%
  mutate(., smpl = variable) %>%
  mutate(., Class1 = Class) %>%
  mutate(., Ion1 = Ion) %>%
  group_by(Class, Ion, variable) %>%
  group_modify(~ {
    .x %>%
      mutate(., IS_int = as.numeric(GPMV.IS.w.conc[GPMV.IS.w.conc$Class == unique(Class1) & GPMV.IS.w.conc$Ion == unique(Ion1), unique(smpl)])) %>%
      mutate(., conc_in_smpl = as.numeric(GPMV.IS.w.conc[GPMV.IS.w.conc$Class == unique(Class1) & GPMV.IS.w.conc$Ion == unique(Ion1), "Conc_in_sample"])) %>%
      mutate(., Endo_conc = value/(IS_int/conc_in_smpl))
  })

df.conc.wide <- df.conc %>%
  select(., c("LipidIon", "Endo_conc")) %>%
  pivot_wider(names_from = variable, values_from = Endo_conc)

GPMV.w.conc <- df.conc.wide
save(GPMV.w.conc, file = paste0("./Output_Data/GPMV_483lipid_concentration.Rdata"))

lipid.no.conc <- lipid.no.conc %>%
  ungroup(FattyAcid) %>%
  select(., -("FattyAcid"))

GPMV.565.lpd.all <- bind_rows(GPMV.w.conc, lipid.no.conc)
save(GPMV.565.lpd.all, file = paste0("./Output_Data/GPMV_All_483_conc_82_no_conc.Rdata"))
