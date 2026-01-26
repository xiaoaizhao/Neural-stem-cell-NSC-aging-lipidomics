
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(ggrepel)

load("./Output_Data/DB.MsD.Age.ES_Exp2_all_KO.Rdata")
load("./Output_Data/DBPct.MsD.Invitro_Qui_Age_ES.Rdata")
load("./Output_Data/DB.2outof2LC.invitro.Padj.Sig.features.Rdata")

Exp2.ctrl.conc <- Exp2.MsD.DB.es.g.allKO %>% 
  filter(KO == "N") 

E2 <- Exp2.ctrl.conc %>% 
  select(-KO) %>% 
  mutate(Exp = "In vitro Experiment #2")

Invitro <- Invitro.Qui.MsD.CONC.DB.es.g %>% 
  mutate(Exp = "In vitro")

LC2 <- bind_rows(E2, Invitro)

LC2$Age_mean_Old <- as.numeric(LC2$Age_mean_Old)
LC2$Age_mean_Young <- as.numeric(LC2$Age_mean_Young)

LC2 <- LC2 %>% 
  rowwise() %>% 
  mutate(DB.avg = rowMeans(across(c(Age_mean_Old, Age_mean_Young))) * 100) %>% 
  mutate(Sig.lbl = ifelse(Cla_DB %in% unique(padj.sig.db.$Cla_DB), "Significant", "Not significant")) %>% 
  mutate(Sig.txt = ifelse(Cla_DB %in% unique(padj.sig.db.$Cla_DB), Cla_DB, ""))

LC2$Sig.lbl <- factor(LC2$Sig.lbl, levels = c("Not significant" , "Significant"))
LC2$Exp <- factor(LC2$Exp, levels = c( "In vitro", "In vitro Experiment #2"))

LC2 <- LC2 %>% 
  arrange(Sig.lbl)

pal <- c("grey80", "orchid4")

a <- ggplot(LC2, aes(reorder(Cla_DB, DB.avg), DB.avg)) 
a+
  geom_point(aes(color = Sig.lbl), alpha = 0.85, size = 3)+
  # geom_point(aes(color = Sig.lbl), alpha = 0.75, size = 3)+
  theme_classic()+
  geom_text_repel(aes(label = Sig.txt), fontface = 'plain',
                  size = 3,colour = "black",
                  box.padding = unit(0.4, "lines"),
                  seed = 1234,
                  min.segment.length = 0,
                  max.overlaps = 55)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(title = "Significant 2 out of 2 LC-MS in vitro DB abundance")+
  scale_x_discrete(expand = c(0.1,0.1)) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Double bond abundance (mol%) in respective class") +
  facet_wrap(~Exp,  nrow = 1)

ggsave(filename = "./Figure_Panels/EDFig.4c.pdf", width = 6, height = 5, useDingbats = FALSE)
