## Need to run the following script first
### 'Fig3.SideChain_Invitro.vs.GPMV.R'

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(ggrepel)
### ==== 39 side chain features are too many to label, add text to only top 30% that go up and down with age ====
### ==== color code all significant SC that are identified in both GPMV and in vitro ====
load("./Output_Data/FDR.sig.SC_GPMV_Invitro.top30.pct.Rdata")
load("./Output_Data/FDR.sig.SC_GPMV_Invitro.Rdata")

load("./Output_Data/Invitro.SC.analysis.Rdata") 

Invitro.SC <- Invitro.SumSC %>% 
  mutate(Cla_SC = paste0(Class, "(", Sep.SC, ")")) %>% 
  mutate(Age = ifelse(grepl("^O", Samples), "Old", "Young")) %>% 
  group_by(Cla_SC) %>% 
  summarise(meanSumSC = mean(SumSC)) %>% 
  mutate(Exp = "In vitro")

load("./Output_Data/GPMV.SC.analysis.Rdata")

GPMV.SC <- GPMV.SC.sum %>% 
  mutate(Cla_SC = paste0(Class, "(", Sep.SC, ")")) %>% 
  mutate(Age = ifelse(grepl("^O", Sample), "Old", "Young")) %>% 
  group_by(Cla_SC) %>% 
  summarise(meanSumSC = mean(SumSC)) %>% 
  mutate(Exp = "GPMV")


GPMV.invitro <- bind_rows(Invitro.SC, GPMV.SC) %>% 
  rowwise() %>% 
  mutate(Sig.lbl = ifelse(Cla_SC %in% SC.sig.GPMV.invitro, "Significant", "Not significant")) %>% 
  mutate(Sig.txt = ifelse(Cla_SC %in% top30.ls, Cla_SC, "")) %>% 
  mutate(Log2.SCSum = log2(meanSumSC)) 

GPMV.invitro$Sig.lbl <- factor(GPMV.invitro$Sig.lbl, levels = c("Not significant" , "Significant"))
GPMV.invitro$Exp <- factor(GPMV.invitro$Exp, levels = c( "GPMV", "In vitro"))

GPMV.invitro <- GPMV.invitro %>% 
  arrange(Sig.lbl)

pal <- c("grey80", "orchid4")

a <- ggplot(GPMV.invitro, aes(reorder(Cla_SC, Log2.SCSum), Log2.SCSum)) 
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
  labs(title = "Significant In vitro and GPMV side chain abundance")+
  scale_x_discrete(expand = c(0.1,0.1)) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Log2 side chain concentration (uM) in respective class") +
  facet_wrap(~Exp,  nrow = 1)

ggsave(filename = "./Figure_Panels/EDFig.8g.pdf", width = 6, height = 5, useDingbats = FALSE)

# t <- LC2 %>% 
#   filter(Sig.lbl == "Significant") %>% 
#   group_by(Exp) %>% 
#   tally()
