### Separate DB PCT and raw lipid data into samples of _Mboat2_ overexpression and young plasma membrane lipid supplementation samples
## --------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(RColorBrewer)
library(rstatix)

setwd(rstudioapi::getActiveProject())
## -------------------------------------------------------------------------------------------------------------------
load('Output_Data/MsD.M2PM_CONC.DB_PCT_by_Class.Rdata') 
load("./Output_Data/M2PM_backtoRAW_MsD_Norm_Impt_397_lipids.Rdata") 

M2OE.DB <- MsD.DB.pct.M2PM %>% 
  filter(grepl("_EGFP|_Mb2_OE", Sample)) #976/16 = 61

PM.sup.DB <- MsD.DB.pct.M2PM %>% 
  filter(grepl("_meth|lpd", Sample)) #1464/24 = 61

save(M2OE.DB, file = "./Output_Data/Mboat2_OE_DB_PCT.Rdata")
save(PM.sup.DB, file = "./Output_Data/PM_supp_DB_PCT.Rdata")

## Separate individual lipid data matrix of 2 studies, Mboat2 OE and GPMV supplementation
## -------------------------------------------------------------------------------------------------------------------

M2OE.Lipid <- M2PM.MsD.lpd.conc.all %>% 
  select(matches("_EGFP|_Mb2_OE")) # 16 samples

PM.sup.Lipid <- M2PM.MsD.lpd.conc.all %>% 
  select(matches("_meth|lpd")) # 24 samples

save(M2OE.Lipid, file = "./Output_Data/Mboat2_OE_LIPID.Rdata")
save(PM.sup.Lipid, file = "./Output_Data/PM_supp_LIPID.Rdata")
