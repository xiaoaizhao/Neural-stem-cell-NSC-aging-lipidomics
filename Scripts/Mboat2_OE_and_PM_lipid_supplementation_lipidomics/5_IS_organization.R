### Processing script 5, _Mboat2_ overexpression and young plasma membrane lipid supplementation lipidomics. Lipid extracts were solubilized in 150μl of 90:10 methanol:toluene
# 1. Organize IS matrix
# 14 IS standards with correct ion, will be used for concentration calculation
setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)

##Input IS concentration 
#New stock Cholesterol (d7) is in a slightly different concentration, 1304.53μM instead of 1336.0600μM
ion.tbl <- read.csv("./Input_Data/Final_ion_list_for_cleanUp_052620_for2017LC.csv", stringsAsFactors = F)
load("./Input_Data/IS_new_conc_July_batch1_2.Rdata")
load("./Output_data/IS_M2PM.Rdata") #from script 0_EquiSplash_combine

## remove endogenous cholesterol data
IS.M2PM.org <- IS.M2PM %>% 
  select(., c("Filename",matches("d7|d9"))) %>%
  rename("Samples" = "Filename") 

IS.M2PM.org.t <- as.data.frame(t(IS.M2PM.org)) 
colnames(IS.M2PM.org.t) <- IS.M2PM.org.t[1,]
IS.M2PM.org.t <- IS.M2PM.org.t[-1,]

IS.M2PM.org.t.df <- IS.M2PM.org.t %>% 
  rownames_to_column(., var = "IS") %>% 
  mutate(., Class = c("LPC", "LPE", "MG","PI", "PS", "PG", "SM", "Cholesterol", "PC", "Cer", "PE", "DG", "TG", "ChE"
  )) %>% 
  relocate(Class, .after = "IS") %>% 
  mutate_at(vars(contains("_")), as.numeric) %>% 
  mutate(., Ion = substr(row.names(IS.M2PM.org.t),
                         stri_locate_last(row.names(IS.M2PM.org.t), regex = "\\+|\\-"), 
                         stri_locate_last(row.names(IS.M2PM.org.t), regex = "_")-1)) %>% 
  mutate(., Ion = ifelse(Class %in% c("Cholesterol","Cer"),
                         "+H-H2O", Ion)) %>%
  mutate(., Ion = paste0("[", Ion, "]")) %>% 
  mutate(IS2 = IS) %>% 
  group_by(IS) %>%
  group_modify(~ {
    .x %>%
      mutate(., ion_match = ifelse(Ion == ion.tbl$Ion[match(Class, ion.tbl$Class)], "T", "F"))
  }) %>% #14
  filter(., ion_match == "T")  %>%
  arrange(., Ion)   #14

IS.M2PM.to.merge <- IS.M2PM.org.t.df %>% 
  ungroup() %>% 
  select(-ion_match) %>% 
  mutate(., ion_sign = substr(Ion, 2, nchar(Ion)-1)) %>% 
  mutate(`Mixture Component` = substr(IS, 1, stri_locate_last(IS, regex = "_")-1)) %>%
  mutate(`Mixture Component` = substr(`Mixture Component`, 1, stri_locate_last(`Mixture Component`, regex = "\\+|\\-")-1)) %>% 
  mutate(`Mixture Component` = ifelse(Class %in% c("Cholesterol","Cer"),
                                      substr(`Mixture Component`, 1, stri_locate_last(`Mixture Component`, regex = "\\+|\\-")-1),
                                      `Mixture Component`)) %>% 
  select(., c(`Mixture Component`, matches("XZ|QC"), Class, ion_sign)) %>% 
  rename(Ion = ion_sign)

IS_w_conc.M2PM <- left_join(IS.M2PM.to.merge, IS.Jul23.Batch1_2, by = "Mixture Component") 

save(IS_w_conc.M2PM, file = paste0("./Output_data/IS_M2PM_w_conc.Rdata"))
