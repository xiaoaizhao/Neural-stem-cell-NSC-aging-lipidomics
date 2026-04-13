### Plot abundance on individual lipids in 2 in vitro studies
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggbreak)

source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Invitro_lpd.mol.pct.Rdata") 
load("./Output_Data/Invitro_Exp2_lpd.mol.pct.Rdata") 

load("./Output_Data/Lipids.2of2LC.invitro.FDRSig.features.final.Rdata")
lpd.invitro.2of2.sig <- unique(lpd.c2.in.2lc.df.sig.t.final$LipidIon) #80
load("./Output_Data/Mboat2.responsive_lipid_list.Rdata")

E2 <- E2.lpd.mol.pct %>% 
  group_by(LipidIon) %>% 
  summarise(Mean.mol.pct = mean(lpd.mol.pct)) %>% 
  mutate(Exp = "In vitro #2")

Invitro <- Invitro.lpd.mol.pct %>% 
  group_by(LipidIon) %>% 
  summarise(Mean.mol.pct = mean(lpd.mol.pct)) %>% 
  mutate(Exp = "In vitro")

all2invitro <- bind_rows(E2, Invitro)


all2 <- all2invitro %>% 
  mutate(Sig.lbl = ifelse(LipidIon %in% lpd.invitro.2of2.sig, "Significant", "Not significant")) %>% 
  mutate(Sig.txt = ifelse(LipidIon %in% res.lpd.ls, LipidIon, ""))

all2$Exp <- factor(all2$Exp, levels = c("In vitro", "In vitro #2"))
all2$Sig.lbl <- factor(all2$Sig.lbl, levels = c("Not significant", "Significant"))

all2 <- all2 %>% 
  arrange(Sig.lbl)

pal <- c("grey80", "grey30")

## ==== Plot in vitro data first ====
Invitro <- all2 %>% 
  filter(Exp == "In vitro")

Invitro$Exp <- as.character(Invitro$Exp)

a <- ggplot(Invitro, aes(reorder(LipidIon, Mean.mol.pct), Mean.mol.pct))
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
  labs(title = paste0(unique(Invitro$Exp), " significant both in vitro"))+
  scale_x_discrete(expand = c(0.1,0.1)) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Lipid mol%") +
  geom_hline(yintercept = median(Invitro$Mean.mol.pct), linetype = "dashed")
ggsave(paste0("./Figure_Panels/EDFig.1f.invitro_1.pdf"), width = 5, height = 5)

## ==== Plot in vitro Exp #2 next ====
exp2 <- all2 %>% 
  filter(Exp == "In vitro #2")

exp2$Exp <- as.character(exp2$Exp)

a <- ggplot(exp2, aes(reorder(LipidIon, Mean.mol.pct), Mean.mol.pct))
a+
  geom_point(aes(color = Sig.lbl), alpha = 0.85, size = 3)+
  scale_y_break(
    breaks    = c(3, 4.3),          # cut from 7 to 17
    scales    = c(1, 0.1),           # 1 = normal, 2 = 2× height (expand upper part)
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
  labs(title = paste0(unique(exp2$Exp), " significant both in vitro"))+
  scale_x_discrete(expand = c(0.1,0.1)) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Lipid mol%") +
  geom_hline(yintercept = median(exp2$Mean.mol.pct), linetype = "dashed")
ggsave(paste0("./Figure_Panels/EDFig.1f.invitro_2.pdf"), width = 5, height = 5)
