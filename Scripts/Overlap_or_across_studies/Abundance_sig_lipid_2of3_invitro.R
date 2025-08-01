##Prerequisite: Need to run `Lipid_aging_invitro_summary_Exp1-Exp3.R` first
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(ggrepel)


source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_Data/Ef_Size_CONC.Lipid_Exp2_all_KO.Rdata") #Exp2
load("./Output_Data/Exp3_Qui_Age_ES.conc.lipids.Rdata") #Exp3


Exp2.ctrl.conc <- Exp2.CONC.lpd.es.g.allKO %>% 
  filter(KO == "N") #612

load("./Output_Data/Lipids.2outof3LC.invitro.Sig.features.Rdata")
lpd.invitro.2of3.sig <- unique(lpd.c2.in.3lc.df.sig$LipidIon) #127

E2 <- ether.rename(Exp2.ctrl.conc) %>% 
  select(-KO) %>% 
  mutate(Exp = "Primary culture #2")

E3 <- E3.Q.conc.AG.ES %>% 
  mutate(Exp = "Primary culture #3")

all3invitro <- bind_rows(E2, E3)

all3invitro$Age_mean_Old <- as.numeric(all3invitro$Age_mean_Old)
all3invitro$Age_mean_Young <- as.numeric(all3invitro$Age_mean_Young)

all3 <- all3invitro %>% 
  rowwise() %>% 
  mutate(Lpd.avg = rowMeans(across(c(Age_mean_Old, Age_mean_Young)))) %>% 
  mutate(Sig.lbl = ifelse(LipidIon %in% lpd.invitro.2of3.sig, "Significant", "Not significant")) %>% 
  mutate(Sig.txt = ifelse(LipidIon %in% lpd.invitro.2of3.sig, LipidIon, ""))

all3$Sig.lbl <- factor(all3$Sig.lbl, levels = c("Not significant", "Significant"))
all3$Exp <- factor(all3$Exp, levels = c("Primary culture #2", "Primary culture #3"))

all3 <- all3 %>% 
  arrange(Sig.lbl)

t <- all3 %>% 
  group_by(Exp, Sig.lbl) %>% 
  tally()

pal <- c("grey80", "orchid4")

all3.p<- all3 %>%
  mutate(., Exp1 = Exp) %>%
  group_by(., Exp1) %>%
  group_map(~ {
    ggplot(.x) + aes(reorder(LipidIon, Lpd.avg), log2(Lpd.avg)) +
    geom_point(aes(color = Sig.lbl), alpha = 0.85, size = 3)+
  theme_classic()+
  geom_text_repel(aes(label = Sig.txt), fontface = 'plain',
                  size = 3,colour = "black",
                  box.padding = unit(0.4, "lines"),
                  seed = 1234,
                  min.segment.length = 0,
                  max.overlaps = 55)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(title = paste0(unique(.x$Exp), "significant 2 of 3 in vitro"))+
  scale_x_discrete(expand = c(0.1,0.1)) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Lipid abundance") +
  geom_hline(yintercept = median(log2(.x$Lpd.avg)), linetype = "dashed")
ggsave(paste0("./Figure_Panels/EDFig.1i.", unique(.x$Exp),".pdf"), width = 5, height = 5)
  }
)

