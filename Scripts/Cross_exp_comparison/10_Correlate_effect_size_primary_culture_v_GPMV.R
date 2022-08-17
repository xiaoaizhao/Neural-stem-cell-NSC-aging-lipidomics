
##Correlation analysis between primary culture #1 (whole cell extract) and GPMV lipidomics data
##Using Hedge's g effect size
##Subset to only include membrane lipid classes and then generate correlation plot.
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(ggrepel)
setwd(rstudioapi::getActiveProject())

load("./Output_Data/Ef_Size_Lipid_Age_qNSC_InVitro_LC.Rdata")
load("./Output_Data/Ef_Size_Lipid_GPMV.Rdata")

Ef.size.GPMV.LC <- inner_join(GPMV.lpd.es.g, InVitro.lpd.es.g, by = "LipidIon") %>%
  select(., matches("LipidIon|es.g")) %>%
  rename(., "Effect_size_GPMV"= "es_g.x") %>%
  rename(., "Effect_size_invitro"= "es_g.y") 

####correlation on membrane lipids only####
clas.list <- c("Cer", "Cholesterol", "Hex1Cer", "LPC", "LPE", "PC", "PE", "PG", "PI", "PS", "SM")
mem.lipids <- Ef.size.GPMV.LC %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  filter(., Class %in% clas.list) #196

####correlation between two datasets, use ggplot with ggrepel, annotate lipids that are in the top or bottom 10% in both datasets
Eff_size_p <- mem.lipids %>% 
  ungroup() %>% 
  mutate(., NoIon = substr(LipidIon, 1, str_locate(LipidIon, "\\)"))) %>% 
  mutate(., TopHit = ifelse(
    Effect_size_GPMV >=quantile(Effect_size_GPMV, 0.9) &
      Effect_size_invitro >=quantile(Effect_size_invitro, 0.9) |
      Effect_size_GPMV <=quantile(Effect_size_GPMV, 0.1) & 
      Effect_size_invitro <=quantile(Effect_size_invitro, 0.1),
    NoIon, ""
  ))


a <- ggplot(Eff_size_p, aes(Effect_size_GPMV, Effect_size_invitro))
a+geom_smooth(method = 'lm', alpha = 0.3, colour = "blue4", fill='lightskyblue',linetype=1)+
  geom_point(size = 2.5, colour = "grey25", alpha = 0.8)+
  stat_cor(method = "pearson")+
  theme_classic()+
  labs(title = "Effect Size Correlation" , 
       x = "GPMV", 
       y = "Primary NSC Culture")+
  geom_text_repel(aes(label = TopHit), fontface = 'plain',
                  size = 3,colour = "black",
                  xlim = c(min(Eff_size_p$Effect_size_GPMV)-0.2, 
                           max(Eff_size_p$Effect_size_GPMV)+0.2),
                  box.padding = unit(0.55, "lines"),
                  seed = 1234,
                  min.segment.length = 0,
                  max.overlaps = 15)+
  theme(text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(legend.position = "none")
ggsave(filename = "./Figure_Panels/Fig_3c.pdf", height = 5, width = 5, useDingbats=FALSE)
