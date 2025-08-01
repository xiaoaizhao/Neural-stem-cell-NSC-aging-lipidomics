## Overlapping lipids between in vitro summary (2 out of 3 in vitro LC-MS/MS data) with in vivo isolated qNSC lipidomics

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(eulerr)
source("./Scripts/Function_scripts/Effect_size_functions.R")

## In vitro summary
load("./Output_Data/Lipids.2outof3LC.invitro.allfeatures.Rdata")
invitro.sum <- p.c2.in.3lc.df %>% 
  select(LipidIon, MeanES) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
    .x %>% 
      summarise(SumES = unique(MeanES))
  })

## In vivo isolated qNSC lipidomics
load("./Output_Data/Ef_Size_CONC.Lipid_InVivo.Rdata")
Invivo.n <- ether.rename(Invivo.CONC.lpd.es.g)

Invitrosum.Invivo <- list(
  Invitro.sum = invitro.sum$LipidIon,
  Invivo = Invivo.n$LipidIon)


## Venn diagram
pdf(file=paste0("./Figure_Panels/EDFig.5b.left.pdf"), onefile=FALSE)
p1 <- euler(Invitrosum.Invivo)
plot(p1, fills = c("#d9e8e5", "#edd9ba"), 
     quantities = TRUE, legend = FALSE)
dev.off()

Overall.invitrosum.invivo <- list(invitro.sum, Invivo.n) %>% 
  reduce(inner_join, by = "LipidIon") %>% 
  select(matches("LipidIon|SumES|es_g")) %>% 
  rename("Effect_size_Invitro_summary" = "SumES") %>% 
  rename("Effect_size_Invivo" = "es_g") %>% 
  mutate(Category = "Overall overlap")

## Separate lipids that goes up or down with age
### Up with age

up.invitro <- invitro.sum %>% 
  filter(SumES > 0) #175

up.invivo <- Invivo.n %>% 
  filter(es_g > 0) #38

Up.old.Invitrosum.Invivo <- list(
  Invitro.sum = up.invitro$LipidIon,
  Invivo = up.invivo$LipidIon)

up.old <- list(up.invitro, up.invivo) %>% 
  reduce(inner_join, by = "LipidIon") %>% 
  select(matches("LipidIon|SumES|es_g")) %>% 
  rename("Effect_size_Invitro_summary" = "SumES") %>% 
  rename("Effect_size_Invivo" = "es_g") %>% 
  mutate(Category = "Overlap - increase with age")
## Venn diagram

pdf(file=paste0("./Figure_Panels/EDFig.5b.middle.pdf"), onefile=FALSE)
p1 <- euler(Up.old.Invitrosum.Invivo)
plot(p1, fills = c("#d9e8e5", "#edd9ba"), 
     quantities = TRUE, legend = FALSE)
dev.off()


### Down with age
dw.invitro <- invitro.sum %>% 
  filter(SumES < 0) #126

dw.invivo <- Invivo.n %>% 
  filter(es_g < 0) #83

Down.old.Invitrosum.Invivo <- list(
  Invitro.sum = dw.invitro$LipidIon,
  Invivo = dw.invivo$LipidIon)


## Venn diagram
pdf(file=paste0("./Figure_Panels/EDFig.5b.right.pdf"), onefile=FALSE)
p1 <- euler(Down.old.Invitrosum.Invivo)
plot(p1, fills = c("#d9e8e5", "#edd9ba"), 
     quantities = TRUE, legend = FALSE)
dev.off()

up.yng <- list(dw.invitro, dw.invivo) %>% 
  reduce(inner_join, by = "LipidIon") %>% 
  select(matches("LipidIon|SumES|es_g")) %>% 
  rename("Effect_size_Invitro_summary" = "SumES") %>% 
  rename("Effect_size_Invivo" = "es_g") %>% 
  mutate(Category = "Overlap - decrease with age")

Invivo.invitro.overlap <- bind_rows(Overall.invitrosum.invivo, up.old, up.yng)
save(Invivo.invitro.overlap, file = "./Output_Data/Invivo.vs.Invitro.summary_overlap.Rdata")
