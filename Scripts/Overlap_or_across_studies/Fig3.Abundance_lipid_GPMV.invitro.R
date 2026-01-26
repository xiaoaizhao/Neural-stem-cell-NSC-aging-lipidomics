
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggbreak)

source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Invitro_lpd.mol.pct.Rdata") 
load("./Output_Data/GPMV_lpd.mol.pct.Rdata") 

load("./Output_Data/Lipids.GPMV.Invitro.FDRSig.features.Rdata")
sig.GPMV.invitro <- unique(GPMV.Invitro.lpd.fdr.sig$LipidIon)

up.ls <- c("LPC(O-16:0)", "LPC(O-18:1)","PI(18:0_20:4)", "PE(P-18:1_18:1)")
down.ls <- c("PI(18:0_20:3)", "PC(16:0_20:3)")
ls.all <- c(up.ls, down.ls)

GPMV <- GPMV.lpd.mol.pct %>% 
  group_by(LipidIon) %>% 
  summarise(Mean.mol.pct = mean(lpd.mol.pct)) %>% 
  mutate(Exp = "GPMV")

Invitro <- Invitro.lpd.mol.pct %>% 
  group_by(LipidIon) %>% 
  summarise(Mean.mol.pct = mean(lpd.mol.pct)) %>% 
  mutate(Exp = "In vitro")

GPMV.invitro <- bind_rows(GPMV, Invitro)


all2 <- GPMV.invitro %>% 
  mutate(Sig.lbl = ifelse(LipidIon %in% sig.GPMV.invitro, "Significant", "Not significant")) %>% 
  mutate(Sig.lbl = ifelse(LipidIon %in% ls.all, LipidIon, Sig.lbl)) %>% 
  mutate(Sig.txt = ifelse(LipidIon %in% ls.all, LipidIon, ""))

all2$Sig.lbl <- factor(all2$Sig.lbl, levels = c("Not significant", "Significant", "LPC(O-16:0)", "LPC(O-18:1)","PI(18:0_20:4)", "PE(P-18:1_18:1)", "PI(18:0_20:3)", "PC(16:0_20:3)"))
all2$Exp <- factor(all2$Exp, levels = c("In vitro", "GPMV"))


all2 <- all2 %>% 
  arrange(Sig.lbl)

pal <- c("grey80", "grey30", "#CC66CC", "#FA2FB7", "#9D349D", "#CC6699", "#00CDCD", "#0067CD")

## ==== Plot in vitro data first ====
Invitro.p <- all2 %>% 
  filter(Exp == "In vitro")

Invitro.p$Exp <- as.character(Invitro.p$Exp)

a <- ggplot(Invitro.p, aes(reorder(LipidIon, Mean.mol.pct), Mean.mol.pct))
a+
  geom_point(aes(color = Sig.lbl), alpha = 0.85, size = 3)+
  scale_y_break(
    breaks    = c(2, 10),          # cut from 7 to 17
    scales    = c(1, 0.05),           # 1 = normal, 2 = 2× height (expand upper part)
    expand    = TRUE,              # or FALSE / expansion(mult = c(0.05, 0.1))
    space     = 0.15               # gap size between panels (in cm)
  )+
  theme_classic() +
  geom_text_repel(aes(label = Sig.txt, color = Sig.lbl), fontface = 'plain',
                  size = 3,
                  #colour = "black",
                  box.padding = unit(0.4, "lines"),
                  seed = 1234,
                  min.segment.length = 0,
                  max.overlaps = 100)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(title = paste0(unique(Invitro.p$Exp), " significant GPMV + in vitro"))+
  scale_x_discrete(expand = c(0.1,0.1)) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Lipid mol%") +
  geom_hline(yintercept = median(Invitro.p$Mean.mol.pct), linetype = "dashed")
ggsave(paste0("./Figure_Panels/EDFig.8e.", unique(Invitro.p$Exp),".pdf"), width = 5, height = 5)

## ==== Plot in vitro GPMV next ====
GPMV.p <- all2 %>% 
  filter(Exp == "GPMV")

GPMV.p$Exp <- as.character(GPMV.p$Exp)

a <- ggplot(GPMV.p, aes(reorder(LipidIon, Mean.mol.pct), Mean.mol.pct))
a+
  geom_point(aes(color = Sig.lbl), alpha = 0.85, size = 3)+
  scale_y_break(
    breaks    = c(2.5, 10),          # cut from 7 to 17
    scales    = c(1, 0.05),           # 1 = normal, 2 = 2× height (expand upper part)
    expand    = TRUE,              # or FALSE / expansion(mult = c(0.05, 0.1))
    space     = 0.15               # gap size between panels (in cm)
  )+
  theme_classic() +
  geom_text_repel(aes(label = Sig.txt, color = Sig.lbl), fontface = 'plain',
                  size = 3,
                  #colour = "black",
                  box.padding = unit(0.4, "lines"),
                  seed = 1234,
                  min.segment.length = 0,
                  max.overlaps = 100)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(title = paste0(unique(GPMV.p$Exp), " significant GPMV + in vitro"))+
  scale_x_discrete(expand = c(0.1,0.1)) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Lipid mol%") +
  geom_hline(yintercept = median(GPMV.p$Mean.mol.pct), linetype = "dashed")
ggsave(paste0("./Figure_Panels/EDFig.8e.", unique(GPMV.p$Exp),".pdf"), width = 5, height = 5)
