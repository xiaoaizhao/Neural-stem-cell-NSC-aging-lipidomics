### Processing script 5, Primary culture #3 data. Lipid extracts were solubilized in **120μl** of 90:10 methanol:toluene
# 1. Organize IS matrix
# 14 IS standards with correct ion, will be used for concentration calculation
setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)

##Input IS concentration 
#New stock Cholesterol (d7) is in a slightly different concentration, 1304.53μM instead of 1336.0600μM. 
#The final re-solubilized volumn is different between runs as well. Batch 3 samles were solubilized in 120μl. Generate a new data matrix to reflect this difference

IS.label <- read.csv("./Input_Data/Labelled_IS_concentration_2020_052020.csv", stringsAsFactors = F,
                     check.names = F, colClasses = c("character", rep("numeric", 4)))

IS.2023.B3<- IS.label %>% 
  mutate(., `Conc. (uM)` = ifelse(`Mixture Component` == "Cholesterol(d7)", 1304.53, `Conc. (uM)`)) %>% 
  mutate(., Conc_in_Batch3 = ifelse(
    `Mixture Component` == "Cholesterol(d7)", 
    `Conc. (uM)` *0.2/120,
    `Conc. (uM)` *1/120
  ))

save(IS.2023.B3, file = "./Output_data/IS_new_conc_Exp3.Rdata")

### Now organize data matrix of all deuterated IS standards
rm(list = ls())
ion.tbl <- read.csv("./Input_Data/Final_ion_list_for_cleanUp_052620_for2017LC.csv", stringsAsFactors = F)
load("./Output_data/IS_new_conc_Exp3.Rdata")
load("./Output_data/IS_Exp_3.Rdata") #from script 0_EquiSplash_combine

## remove endogenous cholesterol data
IS.Exp3.org <- IS.batch3 %>% 
  select(., c("Filename",matches("d7|d9"))) %>%
  rename("Samples" = "Filename") 

IS.Exp3.org.t <- as.data.frame(t(IS.Exp3.org)) 
colnames(IS.Exp3.org.t) <- IS.Exp3.org.t[1,]
IS.Exp3.org.t <- IS.Exp3.org.t[-1,]


IS.Exp3.org.t <- IS.Exp3.org.t %>% 
  rownames_to_column(., var = "IS") %>% 
  mutate(., Class = c("LPC", "LPE", "MG","PI", "PS", "PG", "SM", "Cholesterol", "PC", "Cer", "PE", "DG", "TG", "ChE"
  )) %>% 
  relocate(Class, .after = "IS") %>% 
  mutate_at(vars(contains("_")), as.numeric) %>% 
  mutate(., Ion = substr(row.names(IS.Exp3.org.t),
                         stri_locate_last(row.names(IS.Exp3.org.t), regex = "\\+|\\-"), 
                         stri_locate_last(row.names(IS.Exp3.org.t), regex = "_")-1)) %>% 
  mutate(., Ion = ifelse(Class %in% c("Cholesterol","Cer"),
                         "+H-H2O", Ion)) %>%
  mutate(., Ion = paste0("[", Ion, "]")) %>%
  relocate(Ion, .after = "Class") %>% 
  mutate(IS2 = IS) %>% 
  group_by(IS) %>%
  group_modify(~ {
    .x %>%
      mutate(., ion_match = ifelse(Ion == ion.tbl$Ion[match(Class, ion.tbl$Class)], "T", "F"))
  }) %>% #14
  filter(., ion_match == "T")  %>%
  arrange(., Ion)   #14

IS.Exp3.to.merge <- IS.Exp3.org.t %>% 
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

IS_w_conc.Exp3 <- left_join(IS.Exp3.to.merge, IS.2023.B3, by = "Mixture Component") 

save(IS_w_conc.Exp3, file = paste0("./Output_data/IS_Exp3_w_conc.Rdata"))
