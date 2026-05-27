
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

load("./Output_Data/Mboat2.responsive_lipid_list.Rdata")

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
  mutate(Sig.txt = ifelse(LipidIon %in% res.lpd.ls, LipidIon, ""))

all2$Sig.lbl <- factor(all2$Sig.lbl, levels = c("Not significant", "Significant"))
all2$Exp <- factor(all2$Exp, levels = c("In vitro", "GPMV"))


all2 <- all2 %>% 
  arrange(Sig.lbl)

pal <- c("grey80", "grey30")

## ==== Plot in vitro data first ====
Invitro.p <- all2 %>% 
  filter(Exp == "In vitro")

Invitro.p$Exp <- as.character(Invitro.p$Exp)

a <- ggplot(Invitro.p, aes(reorder(LipidIon, Mean.mol.pct), Mean.mol.pct))
a+
  geom_point(aes(color = Sig.lbl), alpha = 0.85, size = 3)+
  scale_y_break(
    breaks    = c(1.5, 2),          # cut from 7 to 17
    scales    = c(1, 0.05),           # 1 = normal, 2 = 2× height (expand upper part)
    expand    = TRUE,              # or FALSE / expansion(mult = c(0.05, 0.1))
    space     = 0.15               # gap size between panels (in cm)
  )+
  theme_classic() +
  geom_text_repel(aes(label = Sig.txt), fontface = 'plain',
                  size = 3,
                  colour = "black",
                  box.padding = unit(0.4, "lines"),
                  seed = 1234,
                  min.segment.length = 0,
                  max.overlaps = 200)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(title = paste0(unique(Invitro.p$Exp), " significant GPMV + in vitro"))+
  scale_x_discrete(expand = c(0.1,0.1)) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Lipid mol%") +
  geom_hline(yintercept = median(Invitro.p$Mean.mol.pct), linetype = "dashed")
ggsave(paste0("./Figure_Panels/fig.S9E.", unique(Invitro.p$Exp),".pdf"), width = 5, height = 5)

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
  geom_text_repel(aes(label = Sig.txt), fontface = 'plain',
                  size = 3,
                  colour = "black",
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
ggsave(paste0("./Figure_Panels/fig.S9E.", unique(GPMV.p$Exp),".pdf"), width = 5, height = 5)

unique(GPMV.p$LipidIon)[unique(GPMV.p$LipidIon) %in% res.lpd.ls]
# [1] "LPC(O-16:0)"   "LPC(O-18:1)"   "PC(16:0_22:5)" "PC(17:0_18:1)" "PC(17:1_18:1)" "PE(16:0_16:0)" "PE(18:1_18:2)"
# [8] "PE(18:1_20:4)"