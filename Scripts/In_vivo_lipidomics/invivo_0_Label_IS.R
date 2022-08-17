## Very first script before analyzing _In vivo_ sorted NSC lipidomics data
## Manually annotate all the spike-in standards and save as a separate data frame for downstream normalization and analysis
## Spike-in standards were annotated based on 
# 1. exact mass 
# 2. Mode from which specific ion was detected
# 3. and retention time.

rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
all.df <- read.csv("./Input_Data/Xiaoai_Low_LS_th.csv", stringsAsFactors = F)
##Positive mode first####
all.df <- all.df %>%
  mutate(., IS = "")
Pos <- all.df %>%
  filter(., Mode == "pos")

#DG
DG <- Pos %>%
  filter(., m.z < 605.6 & m.z >605.57) %>%
  mutate(., IS = "15:0-18:1(d7)_DAG+NH4")

#PC
PC <- Pos %>%
  filter(., m.z < 753.7 & m.z >753.6) %>%
  mutate(., IS = "15:0-18:1(d7)_PC+H")

#PE
PE <- Pos %>%
  filter(., m.z < 711.57 & m.z >711.56) %>%
  mutate(., IS = "15:0-18:1(d7)_PE+H")

#TG
TG <- Pos %>%
  filter(., m.z < 829.8 & m.z >829.79) %>%
  mutate(., IS = "15:0-18:1(d7)-15:0_TAG+NH4")

#LPC
LPC <- Pos %>%
  filter(., m.z < 529.4 & m.z >529.39) %>%
  mutate(., IS = "18:1(d7)_Lyso_PC+H")

#Cer
Cer <- Pos %>%
  filter(., m.z < 513.6 & m.z >513.4) %>%
  mutate(., IS = "d18:1-15:0(d7)_Cer+H-H2O")

#SM
SM <- Pos %>%
  filter(., m.z < 738.7 & m.z >738.6 & Retention.time..min.>9) %>%
  mutate(., IS = "d18:1-18:1(d9)_SM+H")

#Cholesterol
Chol <- Pos %>%
  filter(., m.z < 376.4 & m.z >376.3 & Retention.time..min.>9) %>%
  mutate(., IS = "Cholesterol(d7)+H-H2O")

#CE
CE <- Pos %>%
  filter(., m.z < 675.7 & m.z >675.6 & Retention.time..min.>7) %>%
  mutate(., IS = "18:1(d7)_Chol_Ester+NH4")


##Negative mode####
Neg <- all.df %>%
  filter(., Mode == "neg")

#PG
PG <- Neg %>%
  filter(., m.z < 740.6 & m.z >740.5 & Retention.time..min.>9 & Retention.time..min.<10) %>%
  mutate(., IS = "15:0-18:1(d7)_PG-H")  

#PI
PI <- Neg %>%
  filter(., m.z < 828.6 & m.z >828.5 & Retention.time..min.<10) %>%
  mutate(., IS = "15:0-18:1(d7)_PI-H")  

#PS
PS <- Neg %>%
  filter(., m.z < 753.6 & m.z >753.5 & Retention.time..min.<10) %>%
  mutate(., IS = "15:0-18:1(d7)_PS-H") 

#LPE
LPE<- Neg %>%
  filter(., m.z < 485.4 & m.z >485.3 & Retention.time..min.<3) %>%
  mutate(., IS = "18:1(d7)_Lyso_PE-H") 

All.IS <- bind_rows(Cer, Chol, CE, DG, LPC, LPE, PC, PE, PG, PI, PS, SM, TG) %>%
  relocate(., IS, .before = LipidIon)

save(All.IS, file = paste0("./Output_Data/IS_all_sample.Rdata"))

