
## Effect size calculation using Hedge's g on In vivo data
## Effect size of individual lipids
rm(list = ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")

#### Effect size on individual lipid species #########################################################################################
##calculate effect size on In vivo data####

load(file = "./Output_Data/Invivo_Norm_Impt_backtoraw_all130_lipids.Rdata")  #this data is not transformed
all.lipid <- Invivo.raw_int.invivo %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") #1560

all.lipid.df <- all.lipid %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))
Invivo.lpd.es.g <- es.g.func(all.lipid.df, LipidIon, Age, Conc_Int, Sample) %>%
  mutate(., LipidIon = str_replace_all(LipidIon, "/", "_"))

save(Invivo.lpd.es.g, file = paste0("./Output_Data/Ef_Size_Lipid_InVivo.Rdata"))

#### Effect size on individual lipid species - only include lipid with quantitative standards ################################################
##calculate effect size on In vivo data####

load(file = "./Output_Data/Invivo_Norm_Impt_log2_conc_121_lipids.Rdata")  #this data is not transformed
conc.lipid <- 2^Invivo.Impt_norm_conc_all %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") #1452

conc.lipid.df <- conc.lipid %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

Invivo.CONC.lpd.es.g <- es.g.func(conc.lipid.df, LipidIon, Age, Conc, Sample) %>%
  mutate(., LipidIon = str_replace_all(LipidIon, "/", "_"))

save(Invivo.CONC.lpd.es.g, file = paste0("./Output_Data/Ef_Size_CONC.Lipid_InVivo.Rdata"))


#### Effect size on double bond composition only include lipid with quantitative standards ##############################################
#load double bond %mol calculated from previous scripts
load("./Output_Data/Invivo_DB_PCT_by_Class.Rdata") 
Invivo.DB.PCT <- Invivo.DB.PCT %>%
  mutate(., Cla_DB = paste0(Class, DB_num))
Invivo.DB.es.g <- es.g.func(Invivo.DB.PCT, Cla_DB, Age, DB_Pct, Sample)
save(Invivo.DB.es.g, file = paste0("./Output_Data/Ef_Size_DB_pct_Invivo.Rdata"))

#### Effect size on double bond composition #########################################################################################
#load double bond %mol calculated from previous scripts
load("./Output_Data/Invivo_CONC_DB_PCT_by_Class.Rdata") 
Invivo.CONC.DB.PCT <- Invivo.CONC.DB.PCT %>%
  mutate(., Cla_DB = paste0(Class, DB_num))
Invivo.CONC.DB.es.g <- es.g.func(Invivo.CONC.DB.PCT, Cla_DB, Age, DB_Pct, Sample)
Invivo.CONC.DB.es.g.w.CI <- es.g.func.wCI(Invivo.CONC.DB.PCT, Cla_DB, Age, DB_Pct, Sample)
save(Invivo.CONC.DB.es.g, file = paste0("./Output_Data/Ef_Size_CONC_DB_pct_Invivo.Rdata"))
