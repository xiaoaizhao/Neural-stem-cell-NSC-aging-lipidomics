
setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(ComplexHeatmap)
library(GetoptLong)
library(stringi)

source("./Scripts/Function_scripts/Effect_size_functions.R")
load(file = "./Output_data/Spike-in_norm_Mednorm_all_373_lipids_back_to_raw_int.Rdata")

##Subset quiescent NSC samples for statistical testing
qui.lipid <- raw_int %>%
  rownames_to_column(., var = "LipidIon")  %>%
  select(., matches("Qui|LipidIon")) %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Intensity") #4476/373=12 samples

all.lipid.df <- qui.lipid %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))


###Statistical test using wilcoxon rank sum test####
stat.all <- wilcox_stat(all.lipid.df, Intensity, LipidIon)
Exp1.stat <- stat.all
save(Exp1.stat, file = paste0("./Output_data/Exp1_Qui_sample_OvY_stats.Rdata"))