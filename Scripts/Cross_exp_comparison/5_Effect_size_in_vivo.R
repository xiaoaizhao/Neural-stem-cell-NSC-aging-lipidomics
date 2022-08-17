
## Effect size calculation using Hedge's g on In vivo data
## Effect size of individual lipids
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")

#### Effect size on individual lipid species #########################################################################################
##calculate effect size on In vivo data####

load(file = "./Output_Data/Invivo_Norm_Impt_backtoraw_all130_lipids.Rdata")  #this data is not transformed
all.lipid <- raw_int.invivo %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") #1560

all.lipid.df <- all.lipid %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))
Invivo.lpd.es.g <- es.g.func(all.lipid.df, LipidIon, Age, Conc_Int, Sample) %>%
  mutate(., LipidIon = str_replace_all(LipidIon, "/", "_"))

save(Invivo.lpd.es.g, file = paste0("./Output_Data/Ef_Size_Lipid_InVivo.Rdata"))

#### Effect size on double bond composition #########################################################################################
#load double bond %mol calculated from previous scripts
load("./Output_Data/Invivo__DB_PCT_by_Class.Rdata") 
Invivo.DB.PCT <- Invivo.DB.PCT %>%
  mutate(., Cla_DB = paste0(Class, DB_num))
Invivo.DB.es.g <- es.g.func(Invivo.DB.PCT, Cla_DB, Age, DB_Pct, Sample)

save(Invivo.DB.es.g, file = paste0("./Output_Data/Ef_Size_DB_pct_Invivo.Rdata"))

