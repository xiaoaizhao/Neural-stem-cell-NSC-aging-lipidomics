### In vitro #2 - Analysis on overlapping sets of lipids that are increased or decreased in each KO compared to control 

rm(list = ls())
library(tidyverse)
library(ggthemes)
library("scales")
library(ggpubr)
library(UpSetR)
setwd(rstudioapi::getActiveProject())

source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Exp2_278_conc_lipid.MSD.format.Rdata")

KO.l <- E2MsD.frmt %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Concentration") %>% 
  mutate(KO = substr(Sample, nchar(Sample), nchar(Sample))) %>% 
  mutate(Cat = ifelse(KO == "N", "Control", "Target"))

KO.l$KO <- factor(KO.l$KO, levels = c("N", "E", "M", "F", "A", "P"))

ko.target <- c("E", "M", "F", "A", "P")
es.ko <- list()
es.ko1 <- list()

for (ko.t in ko.target) {
  es.ko[[ko.t]] <- KO.l %>% 
    filter(KO == "N" | KO == ko.t)
  df.for.es <- bind_rows(es.ko[[ko.t]])
  es.ko1[[ko.t]] <- es.g.func.KO(df.for.es, LipidIon, KO, Concentration, Sample) %>% 
    select(LipidIon, es_g) %>% 
    mutate(KO = ko.t) %>% 
    mutate(Direction = ifelse(es_g>0, "Up in Ctrl", "Up in KO"))
}

all.ko.es <- bind_rows(es.ko1) 
save(all.ko.es, file = "./Output_Data/Exp2_lpd.es.by.KO.Young+Old.Rdata")

KO.E <- es.ko1$E %>% 
  filter(abs(es_g) >= quantile(abs(es_g), probs = 0.8)) %>% 
  group_split(Direction)

KO.M <- es.ko1$M %>% 
  filter(abs(es_g) >= quantile(abs(es_g), probs = 0.8)) %>% 
  group_split(Direction)

KO.F <- es.ko1$F %>% 
  filter(abs(es_g) >= quantile(abs(es_g), probs = 0.8)) %>% 
  group_split(Direction)

KO.A <- es.ko1$A %>% 
  filter(abs(es_g) >= quantile(abs(es_g), probs = 0.8)) %>% 
  group_split(Direction)

KO.P <- es.ko1$P %>% 
  filter(abs(es_g) >= quantile(abs(es_g), probs = 0.8)) %>% 
  group_split(Direction)
  

#organize matrix for UpSet plot between batch 3 data, primary culture #2 and in vivo
Up.in.Ctrl <- list(
  Elovl5 = paste0(KO.E[[1]]$LipidIon, "_upinCtrl"),
  Fads2 = paste0(KO.F[[1]]$LipidIon, "_upinCtrl"),
  Mboat2 = paste0(KO.M[[1]]$LipidIon, "_upinCtrl"),
  Agpat3 = paste0(KO.A[[1]]$LipidIon, "_upinCtrl"),
  Pla2g4e = paste0(KO.P[[1]]$LipidIon, "_upinCtrl"))

pdf(file= "./Figure_Panels/fig.S12D.bottom.up.in.ctrl.pdf", onefile=FALSE)
upset(fromList(Up.in.Ctrl),nsets = 6,  order.by = c("freq"), 
      sets.bar.color = "#4E6172", 
      mb.ratio = c(0.65, 0.35), 
      text.scale = 1, 
      nintersects = 20)
dev.off()

#organize matrix for UpSet plot between batch 3 data, primary culture #2 and in vivo
Up.in.KO <- list(
  Elovl5 = paste0(KO.E[[2]]$LipidIon, "_upinKO"),
  Fads2 = paste0(KO.F[[2]]$LipidIon, "_upinKO"),
  Mboat2 = paste0(KO.M[[2]]$LipidIon, "_upinKO"),
  Agpat3 = paste0(KO.A[[2]]$LipidIon, "_upinKO"),
  Pla2g4e = paste0(KO.P[[2]]$LipidIon, "_upinKO"))

pdf(file= "./Figure_Panels/fig.S12D.top.up.in.ko.pdf", onefile=FALSE)
upset(fromList(Up.in.KO),nsets = 6,  order.by = c("freq"), 
      sets.bar.color = "#4E6172", 
      mb.ratio = c(0.65, 0.35), 
      text.scale = 1, 
      nintersects = 20)
dev.off()

## === The effects of Mboat2 KO on highlighted lipids ====
### old samples only, use a different effect size calculation to compare between control and KO
KO.l.old <- KO.l %>% 
  mutate(Age = ifelse(grepl("^O", Sample), "Old", "Young")) %>% 
  filter(Age == "Old")

KO.l.old$KO <- factor(KO.l.old$KO, levels = c("N", "E", "M", "F", "A", "P"))

ko.target <- c("E", "M", "F", "A", "P")
es.ko <- list()
es.ko1 <- list()

for (ko.t in ko.target) {
  es.ko[[ko.t]] <- KO.l.old %>% 
    filter(KO == "N" | KO == ko.t) %>% 
    rename("Condition" = "Cat") %>% 
    mutate(Condition = ifelse(KO == "N", "Control", "Treatment"))
  df.for.es <- bind_rows(es.ko[[ko.t]])
  es.ko1[[ko.t]] <- es.g.treat.func(df.for.es, LipidIon, Condition, Concentration, Sample) %>% 
    select(LipidIon, es_g, se_g) %>% 
    mutate(KO = ko.t) %>% 
    mutate(Direction = ifelse(es_g>0, "Up in old KO", "Up in old Ctrl"))
}

all.ko.es.old <- bind_rows(es.ko1) 
save(all.ko.es.old, file = "./Output_Data/Exp2_lpd.es.by.KO.in.old.compare.to.CTRL.Rdata")

## === The effects of Mboat2 KO on side chain analysis====
### old samples only, use a different effect size calculation to compare between control and KO
load("./Output_Data/E2.SC.analysis.Rdata")
KO.SC.old <- E2.SumSC %>% 
  mutate(Age = ifelse(grepl("^O", Samples), "Old", "Young")) %>% 
  filter(Age == "Old") %>% 
  mutate(Cla_SC = paste0(Class, "(", Sep.SC, ")")) %>% 
  mutate(KO = substr(Samples, nchar(Samples), nchar(Samples)))

KO.SC.old$KO <- factor(KO.SC.old$KO, levels = c("N", "E", "M", "F", "A", "P"))

ko.target <- c("E", "M", "F", "A", "P")
es.ko <- list()
es.ko1 <- list()

for (ko.t in ko.target) {
  es.ko[[ko.t]] <- KO.SC.old %>%
    filter(KO == "N" | KO == ko.t) %>%
    mutate(Condition = ifelse(KO == "N", "Control", "Treatment"))
  df.for.es <- bind_rows(es.ko[[ko.t]])
  es.ko1[[ko.t]] <- es.g.treat.func(df.for.es, Cla_SC, Condition, SumSC, Samples) %>% 
    select(Cla_SC, es_g, se_g) %>% 
    mutate(KO = ko.t) %>% 
    mutate(Direction = ifelse(es_g>0, "Up in old KO", "Up in old Ctrl"))
}

SC.all.ko.es.old <- bind_rows(es.ko1) 
save(SC.all.ko.es.old, file = "./Output_Data/Exp2_SC.es.by.KO.in.old.compare.to.CTRL.Rdata")
