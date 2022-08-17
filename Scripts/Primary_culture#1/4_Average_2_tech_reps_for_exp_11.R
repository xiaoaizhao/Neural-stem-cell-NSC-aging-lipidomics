##This lipidomic data set contains samples from two experiments, indicated by "8" or "11" in sample names
##Had two technical replicates for quiescent samples in experiment 11, so take average of the two replicates for downstream analysis

rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

load(file = "./Output_data/Spike-in_normed_373_lipids.Rdata")

df.avg <- norm.lipid.df %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Intensity" ) %>%
  mutate(., ExpDat = ifelse(grepl("11",Sample), "11", "08")) %>%
  mutate(., CellType = ifelse(grepl("qui", Sample), "Quiescent", "Activated")) %>%
  mutate(., ID = substr(Sample, 1, 2)) %>%
  mutate(., SmplCond = paste0(ID, "_", ExpDat, CellType, "_", LipidIon)) ##11160

##get mean intensity from 2 technical replicates (reflecting two different harvesting methods (trypsin + scrape)) of each quiescent samples from "11" experiment, 
##since these are technical replicates from the exact same biological sample
##use average intensity of two samples from each culture of 11 quiescent 
##Keep one row for each lipid (remove the other )
set.seed(12345)
df.avg.long <- df.avg %>%
  group_by(SmplCond) %>%
  group_modify(~ {
    .x %>%
      mutate(., Int_w_11_avg = mean(Intensity)) %>%
      sample_n(., 1) %>%
      mutate(., AvgID = paste0(ID, "_", CellType, "_", ExpDat))
  }) #8952, (11190-8952)/373 = 6, so removed 6 duplicate samples

##pivot matrix back
df.wide <- df.avg.long %>%
  ungroup() %>%
  select(., matches("LipidIon|AvgID|Int_w_11_avg")) %>%
  pivot_wider(names_from = AvgID, values_from =  Int_w_11_avg) ##now instead of 30 samples, there are only 24 samples left


##Save 11 qui averaged data matrix
tg_norm_lipid_373_w_11_avg <- df.wide
save(tg_norm_lipid_373_w_11_avg, file = paste0("./Output_data/Spike-in_norm_373_lipid_w_11Qui_avg.Rdata"))
