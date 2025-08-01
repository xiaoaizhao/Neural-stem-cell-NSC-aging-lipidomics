
# This script for lipidomics processing of _Mboat2_ overexpression and young plasma membrane lipid supplementation
# EquiSplash deuterated standards as well as endogenous cholesterol were manually extracted by Frank using progenesis

rm(list=ls())
setwd(rstudioapi::getActiveProject())
# setwd("~/Dropbox/NSC_aging_lipidomics/Scripts/2023_dataset/Batch_2/")
library(openxlsx)
library(tidyverse)
library(stringi)
library(ggpubr)
### Splash standards detected in both positive and negative mode 
data_book1<-loadWorkbook("./Input_Data/230712_Xiaoai_Lipidomics_B2_Processing_Endogenous_chol_EquiSplash_Labled_Complex_Lipids_extra_ions_230710_Endogenous_Chol.xlsx")
data<-list()

#Row number where samples start and end
start<-5 #check spreadsheet for the beginning of samples
end<-102 #check spreadsheet for the position of the very last sample

for (i in sheets(data_book1)){
  # i="15_0-18_1(d7)-15_0_TAG+NH4"
  data[[i]]<-data.frame(read.xlsx("./Input_Data/230712_Xiaoai_Lipidomics_B2_Processing_Endogenous_chol_EquiSplash_Labled_Complex_Lipids_extra_ions_230710_Endogenous_Chol.xlsx",sheet=i,rows=start:end))
  data[[i]] <- data[[i]] %>% 
    select(., c("Filename", "Area")) %>% 
    rownames_to_column(., var = "Sequence") %>% 
    mutate_at("Sequence", as.numeric) %>%
    mutate(., Mode = ifelse(Sequence < 51, "Positive", "Negative")) %>% #check spreadsheet for the position of the modality transition
    mutate(., IS = paste0(i,"_", Mode)) %>% 
    filter(., !grepl("MSMS", Filename))
  
  data[[i]]$Area <- as.double(data[[i]]$Area)
}

dfclean2 <- data
dfclean2 <- dfclean2[names(dfclean2) != "Component"]

dfsplash2 <- dplyr::bind_rows(dfclean2)

dfsplash2.w <- dfsplash2 %>% 
  select(-c(Mode,Sequence)) %>% 
  pivot_wider(., names_from = "IS", values_from = "Area") %>%
  select(where(~!all(is.na(.x))))

IS.M2PM <- dfsplash2.w %>% 
  select(-c(`d18_1-15_0(d7)_Cer+H_Positive`, 
            `18_1(d7)_MAG+H_Positive`,
            `18_1(d7)_Lyso_PE+H_Positive`))
save(IS.M2PM, file = "./Output_data/IS_M2PM.Rdata")


