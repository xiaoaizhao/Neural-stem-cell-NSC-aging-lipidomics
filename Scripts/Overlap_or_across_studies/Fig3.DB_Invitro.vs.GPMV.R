
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(ggrepel)

load("./Output_Data/DBPct.MsD.Invitro_Qui_Age_ES.Rdata")
load("./Output_Data/DBPct.MsD.GPMV_Age_ES.Rdata")

GPMV <- GPMV.DB.AG.ES %>% 
  mutate(Exp = "GPMV") 

E3 <- Exp3.Qui.MsD.CONC.DB.es.g %>% 
  mutate(Exp = "In vitro")

GPMV.cell <- bind_rows(GPMV, E3)

GPMV.cell$Age_mean_Old <- as.numeric(GPMV.cell$Age_mean_Old)
GPMV.cell$Age_mean_Young <- as.numeric(GPMV.cell$Age_mean_Young)

GPMV.cell <- GPMV.cell %>% 
  rowwise() %>% 
  mutate(DB.avg = rowMeans(across(c(Age_mean_Old, Age_mean_Young))) * 100) %>% 
  mutate(Sig.lbl = ifelse(Cla_DB %in% DB.sig.GPMV.invitro, "Significant", "Not significant")) %>% 
  mutate(Sig.txt = ifelse(Cla_DB %in% DB.sig.GPMV.invitro, Cla_DB, ""))

GPMV.cell$Sig.lbl <- factor(GPMV.cell$Sig.lbl, levels = c("Not significant" , "Significant"))
GPMV.cell$Exp <- factor(GPMV.cell$Exp, levels = c( "GPMV", "In vitro"))

GPMV.cell <- GPMV.cell %>% 
  arrange(Sig.lbl)

pal <- c("grey80", "orchid4")

a <- ggplot(GPMV.cell, aes(reorder(Cla_DB, DB.avg), DB.avg)) 
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
  labs(title = "Significant DB in GPMV and in vitro")+
  scale_x_discrete(expand = c(0.1,0.1)) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Double bond abundance (mol%) in respective class") +
  facet_wrap(~Exp,  nrow = 1)

ggsave(filename = "./Figure_Panels/EDFig.8i.pdf", width = 6, height = 5, useDingbats = FALSE)
