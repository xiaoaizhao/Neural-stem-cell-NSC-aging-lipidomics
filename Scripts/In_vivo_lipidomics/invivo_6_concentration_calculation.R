
##Single-point calculation on complex lipid concentration in sorted qNSC samples
##Each sample is normalized using the spike-in of each class in the same sample
##130 ion cleaned, all have the correct number of side chains, remove duplicated lipids
##121 lipids have a corresponding standard so concentration can be calculated
##9 lipids have no standard, they are lipids in AcCa, LPI and HexCer1 class
rm(list=ls())
library(tidyverse)
library(reshape2)
setwd(rstudioapi::getActiveProject())
load(file = "./Output_Data/Invivo_Ion+FA_clean_130_lipids.Rdata") #130 lipids total
load(file = "./Output_Data/Invivo_13_IS_int_w_conc.Rdata") ##IS value from all sample, with corresponding spike in concentration

lipid.df <- FA_unique130 %>%
  ungroup() %>%
  select(., c("Ion","LipidIon.FA","Class", matches("Y|O", ignore.case = F)))

lipid.no.conc <- lipid.df %>%
  filter(., !Class %in% IS_w_conc$Class) #9 lipids have no IS for quantification

lipid.conc <- lipid.df %>%
  filter(., Class %in% IS_w_conc$Class)  #121 lipids have corresponding IS for quantification

lipidmelt <- melt(lipid.conc) ##1452 rows

lipidmelt$variable <- as.character(lipidmelt$variable)

##append IS value of each sample with corresponding spike in concentration
##single point concentration calculation = [detected_lipid_intensity] /([internal_standard_intensity]/[concentration_of_internal_standard])

df.conc <- lipidmelt %>%
  mutate(., smpl = variable) %>%
  mutate(., Class1 = Class) %>%
  mutate(., Ion1 = Ion) %>%
  group_by(Class, Ion, variable) %>%
  group_modify(~ {
    .x %>%
      mutate(., IS_int = as.numeric(IS_w_conc[IS_w_conc$Class == unique(Class1) & IS_w_conc$Ion == unique(Ion1), unique(smpl)])) %>%
      mutate(., conc_in_smpl = as.numeric(IS_w_conc[IS_w_conc$Class == unique(Class1) & IS_w_conc$Ion == unique(Ion1), "Conc_in_sample"])) %>%
      mutate(., Endo_conc = value/(IS_int/conc_in_smpl))
  })

df.conc.wide <- df.conc %>%
  select(., c("LipidIon.FA", "Endo_conc")) %>%
  pivot_wider(names_from = variable, values_from = Endo_conc)

lipid.w.conc <- df.conc.wide
save(lipid.w.conc, file = paste0("./Output_Data/Invivo_cell_121_lipid_conc.Rdata"))

Cell.lipid.all <- bind_rows(lipid.w.conc, lipid.no.conc)#130
save(Cell.lipid.all, file = paste0("./Output_Data/Invivo_Cell_All_121_conc_9_no_conc.Rdata"))

