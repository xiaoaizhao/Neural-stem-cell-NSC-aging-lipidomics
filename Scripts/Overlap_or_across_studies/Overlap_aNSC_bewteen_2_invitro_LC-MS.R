## Overlapping lipids between 2 LC-MS/MS datasets on primary aNSC cultures

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(eulerr)
source("./Scripts/Function_scripts/Effect_size_functions.R")

load("./Output_Data/Exp3_all.lipid.all.cells.Rdata")
load(file = "./Output_data/Spike-in_norm_Mednorm_all_373_lipids_back_to_raw_int.Rdata")

### Exp1 data
aNSC.Exp1 <- raw_int %>% 
  select(contains("Activated")) %>% 
  rownames_to_column(var = "LipidIon")

A.E1.n <- ether.rename(aNSC.Exp1) %>% 
  rename_all(paste0, "_Exp1") %>% 
  rename("LipidIon" = "LipidIon_Exp1")

### Exp3 data
aNSC.E3 <- Exp3.all.lpd.all.cell %>% 
  select(matches("LipidIon|_aNSC-A"))

A.E3.n <- ether.rename(aNSC.E3) %>% 
  rename_all(paste0, "_Exp3") %>% 
  rename("LipidIon" = "LipidIon_Exp3")

E1E3 <- list(
  Exp1 = A.E1.n$LipidIon,
  Exp3 = A.E3.n$LipidIon)

### Venn diagram
pdf(file=paste0("./Figure_Panels/EDFig.1c.pdf"), onefile=FALSE)
p1 <- euler(E1E3)
plot(p1, fills = c("#eef4fb", "#d9e8e5"),
     quantities = TRUE, legend = FALSE)
dev.off()

Overall.aNSC<- list(A.E1.n, A.E3.n) %>% 
  reduce(inner_join, by = "LipidIon")
save(Overall.aNSC, file = "./Output_Data/aNSC_2_invitro_LCMS_overlap.Rdata")