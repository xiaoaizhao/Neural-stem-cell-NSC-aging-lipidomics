### Effect size calculation on lipids and double bond composition
### Mboat2 overexpression and plasma membrane lipid supplementation

setwd(rstudioapi::getActiveProject())
rm(list=ls())
source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")
library(tidyverse)
#########################################################################################

## ====calculate effect size on Mboat2 OE, all lipids ====####
load("Output_Data/Mboat2_OE_LIPID.Rdata")
M2OE <- M2OE.Lipid %>%  
  rownames_to_column(var = "LipidIon") 
M2OE.df <- MsD.lpd.rmv.abc(M2OE)%>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>% 
  mutate(Condition = ifelse(grepl("_EGFP", Sample), "Control", "Treatment"))  %>% 
  mutate(Age = ifelse(grepl("^O", Sample), "Old", "Young"))

M2OE.old <- M2OE.df %>% 
  filter(Age == "Old")

M2OE.young <- M2OE.df %>% 
  filter(Age == "Young")

M2.OE.ES <- es.g.treat.func(M2OE.df, LipidIon, Condition, Conc, Sample) %>% 
  select(LipidIon, es_g, se_g) 
M2.old.OE.ES <- es.g.treat.func(M2OE.old, LipidIon, Condition, Conc, Sample) %>% 
  select(LipidIon, es_g, se_g)  

# In this dataset, positive effect size means higher in Mboat2 OE, negative effect size means higher in control
save(M2.OE.ES, file = "./Output_Data/M2OE_lpd.es.by.OE.Young+Old.Rdata")
save(M2.old.OE.ES, file = "./Output_Data/M2OE_lpd.es.by.OE.in.Old.Rdata")


## ====calculate effect size on Mboat2 OE, side chain ====####

load("Output_Data/Mboat2OE.SC.analysis.Rdata")

M2OE.SC.old <- M2OE.SC.sum %>% 
  mutate(Cla_SC = paste0(Class, "(", Sep.SC, ")")) %>% 
  mutate(Condition = ifelse(grepl("_EGFP", Sample), "Control", "Treatment"))  %>% 
  mutate(Age = ifelse(grepl("^O", Sample), "Old", "Young")) %>% 
  filter(Age == "Old")

M2.OE.old.SC.ES <- es.g.treat.func(M2OE.SC.old, Cla_SC, Condition, SumSC, Sample) %>% 
  select(Cla_SC, es_g, se_g)  

# In this dataset, positive effect size means higher in Mboat2 OE, negative effect size means higher in control

save(M2.OE.old.SC.ES, file = "./Output_Data/M2OE_SC.es.by.OE.in.Old.Rdata")

