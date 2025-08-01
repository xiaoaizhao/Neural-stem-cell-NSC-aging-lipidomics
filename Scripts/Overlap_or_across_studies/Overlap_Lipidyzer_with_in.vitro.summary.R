## Overlapping lipids between in vitro summary (2 out of 3 in vitro LC-MS/MS data) with Lipidyzer lipidomics
## Note: since annotation output from Lipidyzer is different from LC-MS, we converted lipid annotation from both platforms into a list before overlapping, e.g. c("Cer" ,  "d18:1", "16:0" ) 
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
  }) %>% 
  rowwise() %>% 
  mutate(., ID_string = list(str_split(LipidIon, "\\(|_|\\)")[[1]][1:length(str_split(LipidIon, "\\(|_|\\)")[[1]])-1]))

## Lipidyzer lipidomics on qNSCs
load("./Output_Data/Ldz_Ef_Size_Age_qNSC.Rdata")# Lipidyzer
load("./Output_data/Lipidyzer_qNSC_for.LC-MS.ovlp.Rdata")
Ldz.ovlp <- Lipidyzer.Age.es.g %>% 
  filter(Lipid %in% Ldz.Qui.nodup$LipidIon) %>% 
  rowwise() %>% 
  mutate(ID_string = Ldz.Qui.nodup$ID_string[Ldz.Qui.nodup$LipidIon == Lipid])

Invitrosum.Lipidyzer <- list(
  Invitro.sum = invitro.sum$ID_string,
  Lipidyzer = Ldz.ovlp$ID_string)

## Venn diagram
pdf(file=paste0("./Figure_Panels/EDFig.3c.left.pdf"), onefile=FALSE)
p1 <- euler(Invitrosum.Lipidyzer)
plot(p1, fills = c("#d9e8e5", "#f6f4d6"), 
     quantities = TRUE, legend = FALSE)
dev.off()

Overall.invitrosum.Lipidyzer <- list(invitro.sum, Ldz.ovlp) %>% 
  reduce(inner_join, by = "ID_string") %>% 
  select(matches("LipidIon|SumES|es_g")) %>% 
  rename("Effect_size_Invitro_summary" = "SumES") %>% 
  rename("Effect_size_Lipidyzer" = "es_g") %>% 
  mutate(Category = "Overall overlap")

## Separate lipids that goes up or down with age
### Up with age

up.invitro <- invitro.sum %>% 
  filter(SumES > 0) #175

up.lipidyzer <- Ldz.ovlp %>% 
  filter(es_g > 0) #169

Up.old.Invitrosum.Lipidyzer <- list(
  Invitro.sum = up.invitro$ID_string,
  Lipidyzer = up.lipidyzer$ID_string)

## Venn diagram

pdf(file=paste0("./Figure_Panels/EDFig.3c.middle.pdf"), onefile=FALSE)
p1 <- euler(Up.old.Invitrosum.Lipidyzer)
plot(p1, fills = c("#d9e8e5", "#f6f4d6"), 
     quantities = TRUE, legend = FALSE)
dev.off()

up.old <- list(up.invitro, up.lipidyzer) %>% 
  reduce(inner_join, by = "ID_string") %>% 
  select(matches("LipidIon|SumES|es_g")) %>% 
  rename("Effect_size_Invitro_summary" = "SumES") %>% 
  rename("Effect_size_Lipidyzer" = "es_g") %>% 
  mutate(Category = "Overlap - increase with age")

### Down with age
dw.invitro <- invitro.sum %>% 
  filter(SumES < 0) #126

dw.lipidyzer <- Ldz.ovlp %>% 
  filter(es_g < 0) #128

Down.old.Invitrosum.Lipidyzer <- list(
  Invitro.sum = dw.invitro$ID_string,
  Lipidyzer = dw.lipidyzer$ID_string)


## Venn diagram
pdf(file=paste0("./Figure_Panels/EDFig.3c.right.pdf"), onefile=FALSE)
p1 <- euler(Down.old.Invitrosum.Lipidyzer)
plot(p1, fills = c("#d9e8e5", "#f6f4d6"), 
     quantities = TRUE, legend = FALSE)
dev.off()

up.yng <- list(dw.invitro, dw.lipidyzer) %>% 
  reduce(inner_join, by = "ID_string") %>% 
  select(matches("LipidIon|SumES|es_g")) %>% 
  rename("Effect_size_Invitro_summary" = "SumES") %>% 
  rename("Effect_size_Lipidyzer" = "es_g") %>% 
  mutate(Category = "Overlap - decrease with age")

Lipidyzer.invitro.overlap <- bind_rows(Overall.invitrosum.Lipidyzer, up.old, up.yng)
save(Lipidyzer.invitro.overlap, file = "./Output_Data/Lipidyzer.vs.Invitro.summary_overlap.Rdata")