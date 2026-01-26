### Calculate the mol% of each individual lipid 

setwd(rstudioapi::getActiveProject())
rm(list=ls())
source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")
library(tidyverse)
#########################################################################################
load("Output_Data/GPMV_MsD.Raw.conc.Norm_Impt_291_lipids.Rdata")
GPMV.df <- GPMV.impt.norm.raw.conc.lpd.MsD %>% 
  rownames_to_column(var = "LipidIon")

GPMV <- MsD.lpd.rmv.abc(GPMV.df) %>% 
  pivot_longer(-LipidIon, values_to = "Conc", names_to = "Samples")

total.conc.GPMV <- GPMV %>% 
  group_by(Samples) %>% 
  summarise(TotalLipid = sum(Conc))

GPMV.lpd.mol.pct <- left_join(GPMV, total.conc.GPMV, by = "Samples") %>% 
  mutate(lpd.mol.pct = (Conc/TotalLipid)*100 )

save(GPMV.lpd.mol.pct, file = "./Output_Data/GPMV_lpd.mol.pct.Rdata")
