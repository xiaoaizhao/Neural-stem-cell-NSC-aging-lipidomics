##Prerequisite: Need to run `DB_aging_invitro_summary_Exp1-Exp3.R` first
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(ggrepel)

load("./Output_Data/Ef_Size_DB_pct_InVitro.Rdata") # Primary culture #1
load("./Output_Data/Ef_Size_CONC.DB_PCT_Exp2_all_KO.Rdata") # Primary culture #2
load("./Output_Data/Ef_Size_CONC.DB_pct_Exp3_Qui_aging.Rdata") # Primary culture #3
load("./Output_Data/DB.2outof3LC.invitro.Sig.features.Rdata")

Exp2.ctrl.conc <- Exp2.CONC.DB.es.g.allKO %>% 
  filter(KO == "N") 

E2 <- Exp2.ctrl.conc %>% 
  select(-KO) %>% 
  mutate(Exp = "Primary culture #2")

E1 <- LC.Invitro.DB.es.g %>% 
  mutate(Exp = "Primary culture #1")

E3 <- Exp3.Qui.CONC.DB.es.g %>% 
  mutate(Exp = "Primary culture #3")

LC3 <- bind_rows(E1, E2, E3)

LC3$Age_mean_Old <- as.numeric(LC3$Age_mean_Old)
LC3$Age_mean_Young <- as.numeric(LC3$Age_mean_Young)

LC3 <- LC3 %>% 
  rowwise() %>% 
  mutate(DB.avg = rowMeans(across(c(Age_mean_Old, Age_mean_Young))) * 100) %>% 
  mutate(Sig.lbl = ifelse(Cla_DB %in% p.c2.in.3lc.df.sig$Cla_DB, "Significant", "Not significant")) %>% 
  mutate(Sig.txt = ifelse(Cla_DB %in% p.c2.in.3lc.df.sig$Cla_DB, Cla_DB, ""))

LC3$Sig.lbl <- factor(LC3$Sig.lbl, levels = c("Not significant" , "Significant"))
LC3$Exp <- factor(LC3$Exp, levels = c("Primary culture #1", "Primary culture #2", "Primary culture #3"))

LC3 <- LC3 %>% 
  arrange(Sig.lbl)

t = 25
LC3.tally <- LC3 %>% 
  mutate(Thresh = ifelse(DB.avg>=t, "Abundant", "")) %>% 
  group_by(Exp, Thresh) %>% 
  tally()

pal <- c("grey80", "orchid4")

a <- ggplot(LC3, aes(reorder(Cla_DB, DB.avg), DB.avg)) 
a+
  geom_point(aes(color = Sig.lbl, order = ifelse(Sig.lbl == "Significant", 2, 1)), alpha = 0.85, size = 3)+
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
  labs(title = "Significant 2 out of 3 LC-MS in vitro DB abundance")+
  scale_x_discrete(expand = c(0.1,0.1)) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Double bond abundance (mol%) in respective class") +
  facet_wrap(~Exp,  nrow = 1)

ggsave(filename = "./Figure_Panels/EDFig.4c.pdf", width = 6, height = 5, useDingbats = FALSE)
