##Prerequisite: Need to run `Lipids_GPMV.vs.Exp3.effect.size.correlation.R` first
rm(list = ls())
setwd(rstudioapi::getActiveProject())
library(tidyverse)
library(patchwork)
library(ggpubr)
source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_Data/GPMV.Exp3_overlap_top.bottom.15pct.EfSz.Rdata")

####==== Primary Culture #3 data ====####
load(file = "./Output_data/Exp3_qNSC_quant.lipids.Rdata")

E3.df <- ether.rename(Exp3.conc.Q) %>% 
  pivot_longer(-LipidIon, names_to = "Samples", values_to = "Conc") %>% 
  group_by(LipidIon) %>% 
  summarise(mean.conc = mean(Conc)) %>% 
  mutate(Label15 = ifelse(LipidIon %in% GPMV.Exp3.top15$TopHit, "Top", "Else")) %>% 
  mutate(Label15 = ifelse(LipidIon == "LPC(O-16:2)", "Else", Label15)) %>% 
  arrange(Label15) %>% 
  mutate(Label15.text = ifelse(Label15 == "Top", LipidIon, ""))

E3.df$Label15 <- factor(E3.df$Label15, levels = c("Else", "Top"))
E3.df <- E3.df %>% 
  arrange(Label15)


pal <- c("grey80", "orchid4")
a <- ggplot(E3.df, aes(reorder(LipidIon, mean.conc), log2(mean.conc))) +
  geom_point(aes(color = Label15), alpha = 0.75, size = 3)+
  theme_classic()+
  geom_text_repel(aes(label = Label15.text), fontface = 'plain',
                  size = 3,colour = "black",
                  box.padding = unit(0.4, "lines"),
                  seed = 1234,
                  min.segment.length = 0,
                  max.overlaps = 95)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(title = "Shared lipids in Primary culture #3")+
  scale_x_discrete(expand = c(0.1,0.1)) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Log2 mean concentration (uM)")

####==== GPMV data ====####
source("./Scripts/Function_scripts/Effect_size_functions.R")
load(file = "./Output_Data/GPMV_Norm_Impt_log2_conc_483_lipids.Rdata")
conc.gpmv <- 2^GPMV_Impt_norm_conc_all %>% 
  rownames_to_column(var = "LipidIon")

GPMV.df <- ether.rename(conc.gpmv) %>% 
  mutate(LipidID = LipidIon) %>% 
  select(-LipidIon) %>% 
  pivot_longer(-LipidID, names_to = "Samples", values_to = "Conc") %>% 
  group_by(LipidID) %>% 
  summarise(mean.conc = mean(Conc)) %>% 
  mutate(Label15 = ifelse(LipidID %in% GPMV.Exp3.top15$TopHit, "Top", "Else")) %>% 
  mutate(Label15 = ifelse(LipidID == "LPC(O-16:2)", "Else", Label15)) %>% 
  arrange(Label15) %>% 
  mutate(Label15.text = ifelse(Label15 == "Top", LipidID, ""))

GPMV.df$Label15 <- factor(GPMV.df$Label15, levels = c("Else", "Top"))

GPMV.df <- GPMV.df %>% 
  arrange(Label15)

pal <- c("grey80", "orchid4")
b <- ggplot(GPMV.df, aes(reorder(LipidID, mean.conc), log2(mean.conc))) +
  geom_point(aes(color = Label15), alpha = 0.75, size = 3)+
  theme_classic()+
  geom_text_repel(aes(label = Label15.text), fontface = 'plain',
                  size = 3,colour = "black",
                  box.padding = unit(0.4, "lines"),
                  seed = 1234,
                  min.segment.length = 0,
                  max.overlaps = 150)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(title = "Shared lipids in GPMV")+
  scale_x_discrete(expand = c(0.1,0.1)) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Log2 mean concentration (uM)")

a+b
ggsave(filename = "./Figure_Panels/EDFig.8g.pdf", width = 6, height =4, useDingbats = FALSE)