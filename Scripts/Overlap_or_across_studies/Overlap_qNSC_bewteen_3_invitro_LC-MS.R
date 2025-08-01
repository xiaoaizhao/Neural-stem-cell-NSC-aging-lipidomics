## Overlapping lipids between 3 LC-MS/MS datasets on primary qNSC cultures

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(eulerr)
source("./Scripts/Function_scripts/Effect_size_functions.R")
#### ====import effect size calculation on aging differences of 3 experiments ====
#Exp1
load("./Output_Data/Ef_Size_Lipid_Age_qNSC_InVitro_LC.Rdata")
E1.Q.AG.ES <- ether.rename(InVitro.lpd.es.g)
#Exp2
load("./Output_Data/Ef_Size_Lipid_Exp2_all_KO.Rdata")
E2.ctr.AG.ES <- ether.rename(Exp2.lpd.es.g.allKO) %>% 
  filter(KO == "N")
#Exp3
load("./Output_Data/Exp3_Qui_Age_ES.Rdata")

ls.E1<- E1.Q.AG.ES$LipidIon
ls.E2 <- E2.ctr.AG.ES$LipidIon
ls.E3 <- E3.Q.AG.ES$LipidIon

## ====Overall overlap between 3 studies====
E1E2E3 <- list(
  Exp1 = ls.E1,
  Exp2 = ls.E2,
  Exp3 = ls.E3)


pdf(file=paste0("./Figure_Panels/EDFig.1b.left.pdf"), onefile=FALSE)
p1 <- euler(E1E2E3)
plot(p1, fills = c("#eef4fb", "#d8f0f6", "#d9e8e5"),
     quantities = TRUE, legend = FALSE)
dev.off()

Overall.df <- list(E1.Q.AG.ES, E2.ctr.AG.ES, E3.Q.AG.ES) %>% 
  reduce(inner_join, by = "LipidIon") %>% 
  select(matches("LipidIon|es_g")) %>% 
  rename("Effect_size_Exp1" = "es_g.x") %>% 
  rename("Effect_size_Exp2" = "es_g.y") %>% 
  rename("Effect_size_Exp3" = "es_g") %>% 
  mutate(Category = "Overall overlap")

## ====Overlapping lipids between 3 studies - increase with age====
E1.Q.AG.ES <- E1.Q.AG.ES %>% 
  mutate(Direction = ifelse(es_g > 0 , "Old", "Young")) %>% 
  group_by(Direction) %>% 
  group_split()

Exp1.hi.Old <- E1.Q.AG.ES[[1]]
Exp1.hi.Young <- E1.Q.AG.ES[[2]]

E2.ctr.AG.ES <- E2.ctr.AG.ES %>% 
  mutate(Direction = ifelse(es_g > 0 , "Old", "Young")) %>% 
  group_by(Direction) %>% 
  group_split()

Exp2.hi.Old <- E2.ctr.AG.ES[[1]]
Exp2.hi.Young <- E2.ctr.AG.ES[[2]]

E3.Q.AG.ES <- E3.Q.AG.ES %>% 
  mutate(Direction = ifelse(es_g > 0 , "Old", "Young")) %>% 
  group_by(Direction) %>% 
  group_split()

Exp3.hi.Old <- E3.Q.AG.ES[[1]]
Exp3.hi.Young <- E3.Q.AG.ES[[2]]

hi.old.ls.Exp1<- Exp1.hi.Old$LipidIon
hi.old.ls.Exp2 <- Exp2.hi.Old$LipidIon
hi.old.ls.Exp3<- Exp3.hi.Old$LipidIon

hi.old.venn <- list(
  Exp1 = hi.old.ls.Exp1,
  Exp2 = hi.old.ls.Exp2,
  Exp3 = hi.old.ls.Exp3)

pdf(file=paste0("./Figure_Panels/EDFig.1b.middle.pdf"), onefile=FALSE)
p1 <- euler(hi.old.venn)
plot(p1, fills = c("#eef4fb", "#d8f0f6", "#d9e8e5"),
     quantities = TRUE, legend = FALSE)
dev.off()

Overall.hi.old <- list(Exp1.hi.Old, Exp2.hi.Old, Exp3.hi.Old) %>% 
  reduce(inner_join, by = "LipidIon") %>% 
  select(matches("LipidIon|es_g")) %>% 
  rename("Effect_size_Exp1" = "es_g.x") %>% 
  rename("Effect_size_Exp2" = "es_g.y") %>% 
  rename("Effect_size_Exp3" = "es_g") %>% 
  mutate(Category = "Overlap - increase with age")

## ====Overlapping lipids between 3 studies - decrease with age====
hi.yng.ls.Exp1<- Exp1.hi.Young$LipidIon
hi.yng.ls.Exp2 <- Exp2.hi.Young$LipidIon
hi.yng.ls.Exp3<- Exp3.hi.Young$LipidIon

hi.yng.venn <- list(
  Exp1 = hi.yng.ls.Exp1,
  Exp2 = hi.yng.ls.Exp2,
  Exp3 = hi.yng.ls.Exp3)

pdf(file=paste0("./Figure_Panels/EDFig.1b.right.pdf"), onefile=FALSE)
p1 <- euler(hi.yng.venn)
plot(p1, fills = c("#eef4fb", "#d8f0f6", "#d9e8e5"),
     quantities = TRUE, legend = FALSE)
dev.off()

Overall.hi.yng <- list(Exp1.hi.Young, Exp2.hi.Young, Exp3.hi.Young) %>% 
  reduce(inner_join, by = "LipidIon") %>% 
  select(matches("LipidIon|es_g")) %>% 
  rename("Effect_size_Exp1" = "es_g.x") %>% 
  rename("Effect_size_Exp2" = "es_g.y") %>% 
  rename("Effect_size_Exp3" = "es_g") %>% 
  mutate(Category = "Overlap - decrease with age")

Invitro.overlap <- bind_rows(Overall.df, Overall.hi.old, Overall.hi.yng)
save(Invitro.overlap, file = "./Output_Data/qNSC_3_invitro_LCMS_overlap.Rdata")
