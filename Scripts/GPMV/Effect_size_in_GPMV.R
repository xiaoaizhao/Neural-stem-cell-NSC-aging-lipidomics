
## Effect size calculation using Hedge's g on GPMV data
## Effect size of individual lipids
rm(list = ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")

#### Effect size on individual lipid species - all lipids #####################################################################################
##calculate effect size on GPMV data####

load(file = "./Output_Data/GPMV_Norm_Impt_backtoraw_all565_lipids.Rdata") ##load pre-processed matrix with lipid concentration
all.lipid <- raw_int.GPMV %>%  
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") #9040/565=16

all.lipid.df <- all.lipid %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))
GPMV.lpd.es.g <- es.g.func(all.lipid.df, LipidIon, Age, Conc_Int, Sample) %>%
  mutate(., LipidIon = str_replace_all(LipidIon, "/", "_"))

save(GPMV.lpd.es.g, file = paste0("./Output_Data/Ef_Size_Lipid_GPMV.Rdata"))

#### Effect size on individual lipid species - lipid species with quantitative standards only #################################################
##calculate effect size on GPMV data####
load(file = "./Output_Data/GPMV_Norm_Impt_log2_conc_483_lipids.Rdata")
conc.gpmv <- 2^GPMV_Impt_norm_conc_all %>% 
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc")

conc.gpmv.df <- conc.gpmv %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))
GPMV.conc.lpd.es.g <- es.g.func(conc.gpmv.df, LipidIon, Age, Conc, Sample) %>%
  mutate(., LipidIon = str_replace_all(LipidIon, "/", "_"))

save(GPMV.conc.lpd.es.g, file = paste0("./Output_Data/Ef_Size_Conc.Lipid_GPMV.Rdata"))

#### Effect size on double bond composition #########################################################################################
#load double bond %mol calculated from previous scripts
load("./Output_Data/GPMV_DB_PCT_all_samples.Rdata")
GPMV_DB <- GPMV_DB %>%
  mutate(., Cla_DB = paste0(Class, DB_num))
GPMV.DB.es.g <- es.g.func(GPMV_DB, Cla_DB, Age, DB_Pct, Sample)

save(GPMV.DB.es.g, file = paste0("./Output_Data/Ef_Size_DB_pct_GPMV.Rdata"))

#### Effect size on double bond composition - lipid species with quantitative standards only ##############################################
#load double bond %mol calculated from previous scripts
load("./Output_Data/GPMV_CONC.DB_PCT_all_samples.Rdata")
GPMV_Conc.DB <- GPMV_Conc.DB %>%
  mutate(., Cla_DB = paste0(Class, DB_num))
GPMV.CONC.DB.es.g <- es.g.func(GPMV_Conc.DB, Cla_DB, Age, DB_Pct, Sample)
save(GPMV.CONC.DB.es.g, file = paste0("./Output_Data/Ef_Size_CONC.DB_pct_GPMV.Rdata"))

