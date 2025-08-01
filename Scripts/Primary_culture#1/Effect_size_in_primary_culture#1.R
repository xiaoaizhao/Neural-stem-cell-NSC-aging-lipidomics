
## Effect size calculation using Hedge's g on In vitro LC-MS data
## Effect size of individual lipids as well as on double bond composition
rm(list = ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")

#### Effect size on individual lipid species #########################################################################################
##calculate effect size on in vitro data between age in quiescent cells####
load(file = "./Output_data/Spike-in_norm_Mednorm_all_373_lipids_back_to_raw_int.Rdata")
qui.lipid <- raw_int %>%
  rownames_to_column(., var = "LipidIon")  %>%
  select(., matches("Qui|LipidIon")) %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Intensity") #4476/373=12 samples

invitro.all.lipid.df <- qui.lipid %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

InVitro.lpd.es.g <- es.g.func(invitro.all.lipid.df, LipidIon, Age, Intensity, Sample) %>%
  mutate(., LipidIon = str_replace_all(LipidIon, "/", "_"))

save(InVitro.lpd.es.g, file = paste0("./Output_Data/Ef_Size_Lipid_Age_qNSC_InVitro_LC.Rdata"))



#### Effect size on double bond composition #########################################################################################
#load double bond %mol calculated from previous scripts
load("./Output_Data/Qui_NSC_DB_PCT_by_Class_2017LC-MS.Rdata")

#Double bond effect size between young and old in quiescent cells
Qui_LC2017_df <- Qui_LC2017_df %>%
  mutate(., Cla_DB = paste0(Class, DB_num))

LC.Invitro.DB.es.g <- es.g.func(Qui_LC2017_df, Cla_DB, Age, DB_Pct, Sample)

save(LC.Invitro.DB.es.g, file = paste0("./Output_Data/Ef_Size_DB_pct_InVitro.Rdata"))
