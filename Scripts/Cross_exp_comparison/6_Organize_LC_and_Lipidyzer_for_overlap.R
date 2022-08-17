##Organize LC-MS and Lipidyzer data for overlap analysis
##Used data matrix with effect size because later it will be used as a selection criteria to get rid of duplicated lipids

##Steps:
#1. Get class list from both datasets, intersect to find common classes. 
#2. Separate lipid species base on number of side chains  (i.e, 1 or 2). Since TG side chain is not separated well in lipidyzer, only consider total number of carbon and double bonds when comparing two datasets.
#3. Change LC-MS nomenclature to be consistent with Lipidyzer nomenclature e.g. 18:0e -> O-18:0, 16:0p -> P-16:0
#4. Combine all carbon and double bond detected in 3 chains of TG in LC data into 1 side chain, to match lipidyzer data
#5. Save grouped data frames separately for each dataset. These will be used to generate overlapping lipid list in the next script
rm(list=ls())
library(tidyverse)
library(stringi)

setwd(rstudioapi::getActiveProject())
#### Import both sets of data and determine the overlapping classes of lipids#######################################################################
#### Lipidzyer data
load(file = "./Output_Data/Ldz_Ef_Size_Age_qNSC.Rdata")
Ldz.df <- Lipidyzer.Age.es.g %>%
  rename("LipidIon" = "Lipid") %>%
  mutate(., ID = str_split(LipidIon, "\\.")) %>%
  group_by(LipidIon) %>%
  group_modify(~ {
    .x %>%
      mutate(., Class = ID[[1]][1]) %>%
      mutate(., Class = ifelse(grepl("DAG", Class), "DG", Class)) %>%
      mutate(., Class = ifelse(grepl("TAG", Class), "TG", Class)) %>%
      mutate(., Class = ifelse("CE" %in% Class, "ChE", Class)) 
  }) 

####LC-MS data
load(file = "./Output_Data/Ef_Size_Lipid_Age_qNSC_InVitro_LC.Rdata")
LC.df <- InVitro.lpd.es.g %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1))

common.clas <- intersect(LC.df$Class, Ldz.df$Class)

####Note: Because only one side chains on SM was identified in Lipidzyer data, whereas both side chain is identified on LC-MS. 
####Even though this class of lipid is detected in both datasets, I decide to not include it in overlapping analysis.

clas_commmon_w_lc <- common.clas[!grepl("SM", common.clas)]

#### Organize lipid annotation for both datasets##################################################################################################
#### Lipidzyer data
#uniform class name for DG, TG and ChE, reformat side chain, as well as create a list with class name and FA to compare to LC
ldyz_common_cla_w_LC <- Ldz.df %>%
  filter(., Class %in% clas_commmon_w_lc) %>% ##442 lipids, a majority out of 511 lipids
  mutate(., Lipid1 = LipidIon) %>%
  group_by(LipidIon) %>%
  group_modify(~ {
    .x %>%
      mutate(., FA1 = ifelse(length(ID[[1]]) != 7, paste0(ID[[1]][2], ":", ID[[1]][3]),
                             paste0(ID[[1]][2], "-", ID[[1]][3], ":", ID[[1]][4]))) %>%
      mutate(., FA2 = ifelse(length(ID[[1]]) != 7, paste0(ID[[1]][4], ":", ID[[1]][5]),
                             paste0(ID[[1]][5], ":", ID[[1]][6]))) %>%
      mutate(., FA1 = ifelse("TG" %in% Class, paste0(substr(Lipid1, str_locate(Lipid1, "TAG")+3, str_locate(Lipid1, "TAG")+4),
                                                     ":", ID[[1]][2]), FA1))
  }) 

FA1.df <- ldyz_common_cla_w_LC %>%
  filter(., Class %in% c("LPC", "LPE", "TG")) %>% #230
  mutate(., ID_string = list(c(Class, FA1)))
Dyzer.FA1 <- FA1.df
save(Dyzer.FA1, file = paste0("./Output_data/Lipidyzer_FA1_for_overlap.Rdata"))

FA2.df <- ldyz_common_cla_w_LC %>%
  filter(., !Class %in% c("LPC", "LPE", "TG")) %>% #212
  mutate(., ID_string = list(c(Class, FA1, FA2)))
Dyzer.FA2 <- FA2.df
save(Dyzer.FA2, file = paste0("./Output_data/Lipidyzer_FA2_for_overlap.Rdata"))

#### LC-MS data

LC.df.common.class <- LC.df %>%
  filter(., Class %in% clas_commmon_w_lc)  %>% #232 lipids
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) 
LC.FA1 <- LC.df.common.class %>%
  filter(., Class %in% c("LPC", "LPE")) %>% 
  group_by(LipidIon) %>%
  group_modify(~ {
    .x %>%
      mutate(., ID_string = list(c(Class, SideChain)))
  })

LC.TG <- LC.df.common.class %>%
  filter(., Class %in% "TG")  %>%#41
  group_by(LipidIon) %>%
  group_modify(~ {
    .x %>%
      mutate(., C1 = as.numeric(substr(SideChain, 1, str_locate(SideChain, ":")-1))) %>%
      mutate(., C2 = as.numeric(substr(SideChain, str_locate(SideChain, "_")+1, str_locate_all(SideChain, ":")[[1]][2]-1))) %>%
      mutate(., C3 = as.numeric(substr(SideChain, str_locate_all(SideChain, "_")[[1]][2]+1, str_locate_all(SideChain, ":")[[1]][3]-1))) %>%
      mutate(., db1 = as.numeric(substr(SideChain, str_locate(SideChain, ":")+1, str_locate(SideChain, "_")-1))) %>%
      mutate(., db2 = as.numeric(substr(SideChain, str_locate_all(SideChain, ":")[[1]][2]+1, str_locate_all(SideChain, "_")[[1]][2]-1))) %>%
      mutate(., db3 = as.numeric(substr(SideChain, str_locate_all(SideChain, ":")[[1]][3]+1, nchar(SideChain)))) %>%
      mutate(., Reconstruct_ID = paste0(C1, ":", db1, "_", C2, ":", db2, "_", C3, ":", db3)) %>%
      mutate(., check = ifelse(Reconstruct_ID == SideChain, "T", "F")) %>%
      mutate(., Total_C_db = paste0(sum(C1, C2, C3), ":", sum(db1, db2, db3)))
  }) %>%
  mutate(., ID_string = list(c(Class, Total_C_db)))

TG.to.rbind <- LC.TG %>%
  ungroup() %>%
  select(., -matches("C1|C2|C3|db1|db2|db3|Reconstruct_ID|check|Total_C_db"))

LC.FA1.all <- bind_rows(TG.to.rbind, LC.FA1) #59
save(LC.FA1.all, file = paste0("./Output_data/LC-MS_FA1_for_dyzer_overlap.Rdata"))

LC.FA2 <- LC.df.common.class %>%
  filter(., !Class %in% c("LPC", "LPE", "TG")) %>%  #173
  mutate(., FA1= substr(SideChain, 1, str_locate(SideChain, "_")-1)) %>%
  mutate(., FA2= substr(SideChain, str_locate(SideChain, "_")+1, nchar(SideChain))) %>%
  group_by(LipidIon) %>%
  group_modify(~ {
    .x %>%
      mutate(., FA1_mod = ifelse(grepl("e", FA1), paste0("O-", substr(FA1, 1, nchar(FA1)-1)), FA1)) %>%
      mutate(., FA1_mod = ifelse(grepl("p", FA1_mod), paste0("P-", substr(FA1, 1, nchar(FA1)-1)), FA1_mod)) %>%
      mutate(., FA2_mod = ifelse(grepl("e", FA2), paste0("O-", substr(FA2, 1, nchar(FA2)-1)), FA2)) %>%
      mutate(., FA2_mod = ifelse(grepl("p", FA2_mod), paste0("P-", substr(FA2, 1, nchar(FA2)-1)), FA2_mod)) %>%
      mutate(., ID_string = list(c(Class, FA1_mod, FA2_mod)))
  })

LC.FA2.to.rbind <- LC.FA2 %>%
  ungroup() %>%
  select(., -contains("FA", ignore.case = F))

LC.FA2.all <- LC.FA2.to.rbind
save(LC.FA2.all, file = paste0("./Output_data/LC-MS_FA2_for_dyzer_overlap.Rdata"))

