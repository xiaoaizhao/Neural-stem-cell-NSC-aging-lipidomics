
##Correct lipid with multiple ID - 2nd script following clean duplicated lipid from Second dataset on Primary NSC culture #2
##Notes: Tandem MS/MS was ran on each individual sample, as opposed to ran on QC samples. This is different from the rest of the LC-MS experiment.
##As a result there are no lipids with multiple lipid identification (as seen in Exp2_0_Ion_clean.R). However there are a large number of "isomers"
##I define isomers here as lipids with ppm < 10 and RT within 0.3min. Once they are group, I'll then take one from the group with the highest mscore

####Steps:
#1. Group lipids into isomer group based on ppm < 10 and RT within 0.3min.
#     ppm is parts per million - standard metric for relative delta between two molecules, 
#     ppm is calculated by the difference between any two m/z value and then divide by 10^6
#2. Take one lipid from each isomer group with the highest mscore
#3. Combine lipids after removing isomer, with lipids that are not in any isomer group  - this is the final list of lipids
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

load("./Output_Data/Exp2_FA_clean_738_lipids.Rdata")

org.df <- KO.FA_clean_738 %>%
  select(., matches("LipidIon|RT|CalcMz|mScore|Mean_smple_Int")) %>%
  arrange(., CalcMz) %>%
  mutate(., ID = paste0(LipidIon, "_", RT))

org.df$diff.ppm = NA
org.df$diff.RT = NA

####Group lipids into isomer group based on ppm < 10 and RT within 0.3min.####
isomer_df = list()
for (i in (1:nrow(org.df))) {
  for (j in 1:nrow(org.df)) {
    org.df$diff.ppm[j] <- abs(org.df$CalcMz[j] - org.df$CalcMz[i])/org.df$CalcMz[i] * 10^6
    org.df$diff.RT[j] <- abs(org.df$RT[j] - org.df$RT[i])
  }
  isomer_df[[as.character(org.df$CalcMz[i])]] <- org.df %>%
    filter(., diff.ppm<10 & diff.RT <0.3) %>%
    mutate(., isomer = org.df$CalcMz[i]) #509 isomer groups
}

all.isomer <- bind_rows(isomer_df) #553 lipids included in all isomer groups

non.dup <- org.df %>%
  filter(., !ID %in% all.isomer$ID) #185 lipids does not belong to any isomer group, which are disqualified by ppm, RT or both

####Take one lipid from each isomer group with the highest mscore####
remove.isomer <- all.isomer %>%
  group_by(isomer) %>%
  group_modify(~ {
    .x %>%
      filter(., mScore == max(mScore))
  }) #509 lipids, one for each isomer group

df.rmv.isomer.list <- bind_rows(remove.isomer, non.dup) #694, got rid of 44 lipids 738-694

df.rmv.isomer <- KO.FA_clean_738 %>%
  filter(., Compound %in% df.rmv.isomer.list$ID) #694 after clean up

save(df.rmv.isomer, file = paste0("./Output_Data/Exp2_Dup_ID_rmv+dup_isomer_rmv_694_lipids.Rdata"))
