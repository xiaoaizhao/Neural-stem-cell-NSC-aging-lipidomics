##Input normalization with TG(d5) Spike-in
##Remove 4 Activated NSC pilot sample, they are not unique biogical samples, therefore not part of the study
##30 samples total for this experiment
##Apply normalization factor to lipid data matrix (previous cleaned and filtered, 373 lipids total)
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

TG <- read.csv("./Input_data/TG_(d5)_2017_forR.csv", stringsAsFactors = F)
smpl.list <- read.csv("./Input_data/XZ_lipid_sample_list_031017.csv", stringsAsFactors = F)

TG.smpl <- TG %>%
  select(., matches("Filename|Area")) %>%
  select(., -matches("ISTD.Area")) %>%
  filter(., !grepl("blk|QC", Filename)) %>%
  filter(., grepl("XZ", Filename)) %>%
  mutate(., Smpl_ID = smpl.list$Sample[match(Filename,smpl.list$Key)]) %>%
  filter(., !grepl("Act", Smpl_ID, ignore.case = F)) ##remove 4 Activated NSC pilot sample, they are not part of the study

TG.smpl$Area <- as.numeric(TG.smpl$Area)

##create normalization factor using the median of TG(d5) spike in intensity
TG.IS.norm <- TG.smpl %>%
  mutate(., norm_fctr = as.numeric(scale(TG.smpl$Area, center = F, scale = median(TG.smpl$Area, na.rm = T)))) %>%
  mutate(., Smpl_ID = gsub("\\-", "\\.", Smpl_ID)) ##replace "-" with ".", so that it matches lipid data frame to normalize


##upload cleaned-up dataframe with 373 lipids
load(file = "./Output_data/2_Ion+FA_clean_373_lipids.Rdata")
lipid.df <- FA_unique373 %>%
  ungroup() %>%
  select(., matches("Y|O|LipidIon", ignore.case = F)) %>%
  column_to_rownames(., var = "LipidIon")

#Organize normalization sample list so that it matches the data frame
TG.IS.norm <- TG.IS.norm %>%
  arrange(match(Smpl_ID, colnames(lipid.df)))
all.equal(colnames(lipid.df), TG.IS.norm$Smpl_ID) #[1] TRUE

#Input normalization index
norm.ind <- as.numeric(TG.IS.norm$norm_fctr)

#apply to get normalized data matrix
norm.lipid.df <-sweep(lipid.df, 2, norm.ind, "/")


save(norm.lipid.df, file = paste0("./Output_data/Spike-in_normed_373_lipids.Rdata"))

