
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(eulerr)

source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_Data/Lipids.2outof3LC.invitro.allfeatures.Rdata")
c2.3.all <- p.c2.in.3lc.df %>% 
  select(LipidIon, MeanES) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
    .x %>% 
      summarise(SumES = unique(MeanES))
  })


load("./Output_Data/Ef_Size_Conc.Lipid_GPMV.Rdata")
GPMV.n <- ether.rename(GPMV.conc.lpd.es.g)

Invitrosum.GPMV <- list(
  Invitro.sum = c2.3.all$LipidIon,
  GPMV = GPMV.n$LipidIon)

pdf(file= "./Figure_Panels/EDFig.8c.left.pdf", onefile=FALSE)
p1 <- euler(Invitrosum.GPMV)
plot(p1, fills = c("#d9e8e5", "#E5E5E5"), 
     quantities = TRUE, legend = FALSE)
dev.off()

up.invitro <- c2.3.all %>% 
  filter(SumES > 0) #175

up.gpmv <- GPMV.n %>% 
  filter(es_g > 0) #38

Up.old.Invitrosum.GPMV <- list(
  Invitro.sum = up.invitro$LipidIon,
  GPMV = up.gpmv$LipidIon)

pdf(file= "./Figure_Panels/EDFig.8c.middle.pdf", onefile=FALSE)
p1 <- euler(Up.old.Invitrosum.GPMV)
plot(p1, fills = c("#d9e8e5", "#E5E5E5"), 
     quantities = TRUE, legend = FALSE)
dev.off()

dw.invitro <- c2.3.all %>% 
  filter(SumES < 0) #126

dw.gpmv <- GPMV.n %>% 
  filter(es_g < 0) #83

Down.old.Invitrosum.GPMV <- list(
  Invitro.sum = dw.invitro$LipidIon,
  GPMV = dw.gpmv$LipidIon)

pdf(file= "./Figure_Panels/EDFig.8c.right.pdf", onefile=FALSE)
p1 <- euler(Down.old.Invitrosum.GPMV)
plot(p1, fills = c("#d9e8e5", "#E5E5E5"), 
     quantities = TRUE, legend = FALSE)
dev.off()
