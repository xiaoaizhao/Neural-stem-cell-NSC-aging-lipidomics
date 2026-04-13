### Side chain analysis - summarize the total concentration of lipids of each class that contain the identical side chains
### For SM and Cer, sphingoid base and FA is calculated separately
### For ether lipid, ether-linked side chain is calculated separately from ester-linked side chain

rm(list = ls())
library(tidyverse)
source("./Scripts/Function_scripts/Pre-processing_functions.R")

load("./Output_Data/Exp2_MsD.Norm_Impt_log2_conc_278_lipids.Rdata")

Exp2.raw <- 2^Exp2.log2.impt.norm.conc.MsD %>% 
  rownames_to_column(var = "LipidIon")

E2.rename <- MsD.lpd.format(Exp2.raw) %>% 
  pivot_longer(-LipidIon, values_to = "Conc", names_to = "Samples") %>% 
  rowwise() %>% 
  mutate(SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1)) %>% 
  mutate(Sep.SideChain = ifelse(grepl("^Cer|^SM", LipidIon),
         str_split(SideChain, "/"),
         str_split(SideChain, "_")))

E2.SumSC <- E2.rename %>% 
  mutate(Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>% 
  unnest(Sep.SideChain) %>%
  mutate(Sep.SC = str_trim(Sep.SideChain)) %>%
  group_by(Class, Samples, Sep.SC) %>% 
  summarise(SumSC = sum(Conc))

save(E2.SumSC, file = "./Output_Data/E2.SC.analysis.Rdata")
