
##Quantify classes with standards to get concentration based on one spike-in per class

##Steps:
#1. Make sure that the spike-in standard (IS) used for each class is based on intensity of the same ion as the endogenous lipids
#2. Separate lipids belongs to classes with IS or not
#3. Calculate lipid concentration in classes of lipids with IS - 613 lipids (from 12 classes) 
#   Lipid concentration = [detected_lipid_intensity] /([internal_standard_intensity]/[concentration_of_internal_standard])
#4. Combine with the 81 lipids with intensity only, without concentration
rm(list=ls())
library(tidyverse)
library(reshape2)
setwd(rstudioapi::getActiveProject())

load("./Output_Data/Exp2_Dup_ID_rmv+dup_isomer_rmv_694_lipids.Rdata") #694 lipids, dup iD removed first, then dup isomer removed
load("./Output_Data/Exp2_KO_IS_organized.Rdata")
ion.tbl <- read.csv("./Input_Data/Final_ion_list_for_cleanUp_052620_for2017LC.csv", stringsAsFactors = F)

####Make sure that the spike-in standard used for each class is based on intensity of the same ion as the endogenous lipids####
##MG is removed since the spike-in was detected in a different ion
cla.to.quant <- KO.IS.all %>%
  mutate(., Class1 = Class) %>%
  group_by(Class) %>%
  group_modify(~ {
    .x %>%
      filter(., Ion == as.character(ion.tbl[match(Class1, ion.tbl$Class), "Ion"]))
  }) #12 classes total match determined ion table

####Separate cleaned lipids into lipids with concentration and without concentration ####
lipid.conc <- df.rmv.isomer %>%
  ungroup() %>%
  filter(., Class %in% cla.to.quant$Class1) %>% #613 lipids
  select(., matches("LipidIon|Class|Ion|Y|O", ignore.case = F)) %>%
  select(., -"IonFormula")

##fix column name for classes with spike-in standard
fixcolname <- function(x, na.rm = FALSE){substr(x, nchar(x)-3, nchar(x))}
lipid.conc <- lipid.conc %>%
  rename_at(vars(contains("_")), fixcolname)

##fix column name for classes with spike-in standard
lipid.no.conc <- df.rmv.isomer %>%
  ungroup() %>%
  filter(., !Class %in% cla.to.quant$Class1) %>% #81 lipids
  select(., matches("LipidIon|Class|Ion|Y|O", ignore.case = F)) %>%
  select(., -"IonFormula")
lipid.no.conc <- lipid.no.conc %>%
  rename_at(vars(contains("_")), fixcolname)

####Calculate lipid concentration in classes of lipids with IS####
lipid.melt <- melt(lipid.conc)
lipid.melt$variable <- as.character(lipid.melt$variable)

##single point concentration calculation = [detected_lipid_intensity] /([internal_standard_intensity]/[concentration_of_internal_standard])
df.conc <- lipid.melt %>%
  mutate(., smpl = variable) %>%
  mutate(., Class1 = Class) %>%
  mutate(., Ion1 = Ion) %>%
  group_by(Class, Ion, variable) %>%
  group_modify(~ {
    .x %>%
      mutate(., IS_int = as.numeric(KO.IS.all[KO.IS.all$Class == unique(Class1) & KO.IS.all$Ion == unique(Ion1), unique(smpl)])) %>%
      mutate(., conc_in_smpl = as.numeric(KO.IS.all[KO.IS.all$Class == unique(Class1) & KO.IS.all$Ion == unique(Ion1), "Conc_in_sample"])) %>%
      mutate(., Endo_conc = value/(IS_int/conc_in_smpl))
  })

df.conc.wide <- df.conc %>%
  select(., c("LipidIon", "Endo_conc")) %>%
  pivot_wider(names_from = variable, values_from = Endo_conc)

exp2.lipid.w.conc <- df.conc.wide
save(exp2.lipid.w.conc, file = paste0("./Output_Data/Exp2_613_lipids_w_conc.Rdata"))

exp2.lipid.all <- bind_rows(exp2.lipid.w.conc, lipid.no.conc)
save(exp2.lipid.all, file = paste0("./Output_Data/Exp2_All_613_conc_81_no_conc.Rdata"))


