
# This script is for Primary culture #3 processing
# EquiSplash deuterated standards as well as endogenous cholesterol were manually extracted using Progenesis

rm(list=ls())
setwd(rstudioapi::getActiveProject())
# setwd("~/Dropbox/NSC_aging_lipidomics/Scripts/2023_dataset/Batch_2/")
library(openxlsx)
library(tidyverse)
library(stringi)
library(ggpubr)
### Splash standards detected in both positive and negative mode 
data_book1<-loadWorkbook("./Input_Data/batch3_processing_EquiSplash_Labled_Complex_Lipids_extra_ions_231114_Endogenous_Chol_Chol23.xlsx")
data<-list()

#Row number where samples start and end
start<-5 #check spreadsheet for the beginning of samples
end<-59 #check spreadsheet for the position of the very last sample

for (i in sheets(data_book1)){
  # i="15_0-18_1(d7)-15_0_TAG+NH4"
  data[[i]]<-data.frame(read.xlsx("./Input_Data/batch3_processing_EquiSplash_Labled_Complex_Lipids_extra_ions_231114_Endogenous_Chol_Chol23.xlsx",sheet=i,rows=start:end))
  data[[i]] <- data[[i]] %>% 
    select(., c("Filename", "Area")) %>% 
    rownames_to_column(., var = "Sequence") %>% 
    mutate_at("Sequence", as.numeric) %>%
    mutate(., Mode = ifelse(Sequence < 28, "Positive", "Negative")) %>% #check spreadsheet for the position of the modality transition
    mutate(., IS = paste0(i,"_", Mode)) %>% 
    filter(., !grepl("MSMS", Filename))
  
  data[[i]]$Area <- as.double(data[[i]]$Area)
}

dfclean3 <- data
dfclean3 <- dfclean3[names(dfclean3) != "Component"]

dfsplash3 <- dplyr::bind_rows(dfclean3)

dfsplash3.w <- dfsplash3 %>% 
  select(-c(Mode,Sequence)) %>% 
  pivot_wider(., names_from = "IS", values_from = "Area") %>%
  select(where(~!all(is.na(.x))))

IS.batch3 <- dfsplash3.w %>% 
  select(-c(`d18_1-15_0(d7)_Cer+H_Positive`, 
            `18_1(d7)_MAG+H_Positive`,
            `18_1(d7)_Lyso_PE+H_Positive`))
save(IS.batch3, file = "./Output_data/IS_Exp_3.Rdata")

