### Side chain analysis - summarize the total concentration of lipids of each class that contain the identical side chains
### For SM and Cer, sphingoid base and FA is calculated separately
### For ether lipid, ether-linked side chain is calculated separately from ester-linked side chain

rm(list = ls())
library(tidyverse)
source("./Scripts/Function_scripts/Pre-processing_functions.R")

load("./Output_Data/Invitro_conc.lipid.qNSC.Rdata")

Invitro <- MsD.lpd.format(Invitro.q.conc) %>% 
  pivot_longer(-LipidIon, values_to = "Conc", names_to = "Samples") %>% 
  rowwise() %>% 
  mutate(SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) %>% 
  mutate(Sep.SideChain = ifelse(grepl("^Cer|^SM", LipidIon),
         str_split(SideChain, "/"),
         str_split(SideChain, "_")))

Invitro.SumSC <- Invitro %>% 
  mutate(Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>% 
  unnest(Sep.SideChain) %>%
  mutate(Sep.SC = str_trim(Sep.SideChain)) %>%
  group_by(Class, Samples, Sep.SC) %>% 
  summarise(SumSC = sum(Conc))

save(Invitro.SumSC, file = "./Output_Data/Invitro.SC.analysis.Rdata")
