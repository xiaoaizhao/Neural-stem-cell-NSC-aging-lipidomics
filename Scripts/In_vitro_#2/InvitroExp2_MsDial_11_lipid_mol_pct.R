### Calculate the mol% of each individual lipid 

setwd(rstudioapi::getActiveProject())
rm(list=ls())
source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")
library(tidyverse)
#########################################################################################
load("Output_Data/Exp2_278_conc_lipid.MSD.format.Rdata")

E2.Ctrl <- MsD.lpd.rmv.abc(E2MsD.frmt) %>% 
  select(LipidIon, contains("_N")) %>% 
  pivot_longer(-LipidIon, values_to = "Conc", names_to = "Samples")

total.conc.E2 <- E2.Ctrl %>% 
  group_by(Samples) %>% 
  summarise(TotalLipid = sum(Conc))

E2.lpd.mol.pct <- left_join(E2.Ctrl, total.conc.E2, by = "Samples") %>% 
  mutate(lpd.mol.pct = (Conc/TotalLipid)*100 )

save(E2.lpd.mol.pct, file = "./Output_Data/Invitro_Exp2_lpd.mol.pct.Rdata")
