### Calculate the mol% of each individual lipid 

setwd(rstudioapi::getActiveProject())
rm(list=ls())
source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")
library(tidyverse)
#########################################################################################
load("Output_Data/Invitro_conc.lipid.qNSC.Rdata")

Invitro.Q <- MsD.lpd.rmv.abc(Invitro.q.conc) %>% 
  pivot_longer(-LipidIon, values_to = "Conc", names_to = "Samples")

total.conc.Invitro.Q <- Invitro.Q %>% 
  group_by(Samples) %>% 
  summarise(TotalLipid = sum(Conc))

Invitro.lpd.mol.pct <- left_join(Invitro.Q, total.conc.Invitro.Q, by = "Samples") %>% 
  mutate(lpd.mol.pct = (Conc/TotalLipid)*100 )

save(Invitro.lpd.mol.pct, file = "./Output_Data/Invitro_lpd.mol.pct.Rdata")
