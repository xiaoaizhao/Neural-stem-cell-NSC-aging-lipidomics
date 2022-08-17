
## Correlation analysis between primary NSC culture (2017 LC-MS) data and _In vivo_ data in 2020
## Effect size is calculated as Hedge's g, specifically comparing at lipid species level the effect between age in qNSCs.

## Steps:
#1. Load in effect size matrix of primary culture and data from in vivo sorted cells.
#2. Subset to only keep commonly detected lipids (44 lipids total)
#3. Make correlation plot on effect size between two data sets

##Results: Significant correlation between primary culture and in vivo sorted data set!
rm(list=ls())

library(tidyverse)
library(ggpubr)
library(ggrepel)
setwd(rstudioapi::getActiveProject())

load("./Output_Data/Ef_Size_Lipid_Age_qNSC_InVitro_LC.Rdata")
load("./Output_Data/Ef_Size_Lipid_InVivo.Rdata")

Eff_size_all <- inner_join(Invivo.lpd.es.g, InVitro.lpd.es.g, by = "LipidIon") %>%
  select(., matches("LipidIon|es.g")) %>%
  rename(., "Effect_size_invivo"= "es_g.x") %>%
  rename(., "Effect_size_invitro"= "es_g.y") 



####correlation between two datasets, use ggplot with ggrepel, annotate lipids that are in the top or bottom 10% in either datasets
clas.list <- c("Cer", "Cholesterol", "Hex1Cer", "LPC", "LPE", "PC", "PE", "PG", "PI", "PS", "SM")

Eff_size_all.l <- Eff_size_all %>% 
  ungroup() %>% 
  mutate(., NoIon = substr(LipidIon, 1, str_locate(LipidIon, "\\)"))) %>% 
  mutate(., TopHit = ifelse(
    Effect_size_invivo >=quantile(Effect_size_invivo, 0.8) |
      Effect_size_invitro >=quantile(Effect_size_invitro, 0.8) |
      Effect_size_invivo <=quantile(Effect_size_invivo, 0.2) | 
      Effect_size_invitro <=quantile(Effect_size_invitro, 0.2),
    NoIon, ""
  )) %>% 
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>% 
  mutate(., Lb.mem = ifelse(Class %in% clas.list, TopHit, ""))


a <- ggplot(Eff_size_all.l, aes(Effect_size_invivo, Effect_size_invitro))
a+geom_smooth(method = 'lm', alpha = 0.3, colour = "blue4", fill='lightskyblue',linetype=1)+
  geom_point(size = 3, colour = "grey25")+
  stat_cor(method = "pearson")+
  theme_classic()+
  labs(title = "Effect Size Correlation" , 
    x = "In vivo NSC", 
    y = "Primary NSC Culture")+
  geom_text_repel(aes(label = Lb.mem), fontface = 'plain',
                  size = 3,colour = "black",
                  xlim = c(min(Eff_size_all.l$Effect_size_invivo)-0.2, 
                           max(Eff_size_all.l$Effect_size_invivo)+0.2),
                  box.padding = unit(0.55, "lines"),
                  seed = 1234,
                  min.segment.length = 0,
                  max.overlaps = 10)+
  theme(text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(legend.position = "none")
ggsave(filename = "./Figure_Panels/Fig_1g.pdf", height = 5, width = 5, useDingbats=FALSE)
