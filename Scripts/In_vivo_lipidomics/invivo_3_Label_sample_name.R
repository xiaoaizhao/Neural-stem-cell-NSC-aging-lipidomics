
##change column names to samples names in data matrix
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
load(file = "./Output_Data/Invivo_IonClean_final_159_lipids.Rdata")

label <- read.csv("./Input_Data/Lable_sample_ID.csv", stringsAsFactors = F)

rename.clean <- Ion.clean.159 %>% 
  rename_at(vars(label$Running_ID[1:12]), ~label$Sorted.cell.samples[1:12])

save(rename.clean, file = paste0("./Output_Data/Rename_sample_ion_clean_159_lipids.Rdata"))
