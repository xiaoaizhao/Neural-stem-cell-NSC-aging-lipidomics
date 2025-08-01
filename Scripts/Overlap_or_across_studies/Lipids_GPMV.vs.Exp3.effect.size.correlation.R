
##Correlation analysis between primary culture #3 (whole cell extract) and GPMV lipidomics data
##Using Hedge's g effect size

####====Subset to only include lipids with quantitative standards====####
# 117 overlapping lipids
setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(ggrepel)

load("./Output_Data/Ef_Size_Conc.Lipid_GPMV.Rdata")
load("./Output_Data/Exp3_Qui_Age_ES.conc.lipids.Rdata")

source("./Scripts/Function_scripts/Effect_size_functions.R")
GPMV.es.conc <- ether.rename(GPMV.conc.lpd.es.g)

Ef.size.GPMV.Exp3 <- inner_join(GPMV.es.conc, E3.Q.conc.AG.ES, by = "LipidIon") %>%
  select(., matches("LipidIon|es.g")) %>%
  rename(., "Effect_size_GPMV"= "es_g.x") %>%
  rename(., "Effect_size_invitro"= "es_g.y") 

####correlation between two datasets, use ggplot with ggrepel, annotate lipids that are in the top or bottom 15% in both datasets
Eff_size_p <- Ef.size.GPMV.Exp3 %>%
  ungroup() %>%
  mutate(., TopHit = ifelse(
    Effect_size_GPMV >=quantile(Effect_size_GPMV, 0.85) &
      Effect_size_invitro >=quantile(Effect_size_invitro, 0.85) |
      Effect_size_GPMV <=quantile(Effect_size_GPMV, 0.15) &
      Effect_size_invitro <=quantile(Effect_size_invitro, 0.15),
    LipidIon, ""
  ))

GPMV.Exp3.top15 <- Eff_size_p %>% 
  filter(!TopHit=="")

save(GPMV.Exp3.top15, file = "./Output_Data/GPMV.Exp3_overlap_top.bottom.15pct.EfSz.Rdata")

### after manual validation, among the 9 total 15% lipids, 8 were validated, LPC(O-16:2) was not validated####
# remove LPC(O-16:2)
Eff_size_p.val <- Eff_size_p %>% 
  mutate(TopHit = ifelse(LipidIon == "LPC(O-16:2)", "", TopHit))

a <- ggplot(Eff_size_p.val, aes(Effect_size_GPMV, Effect_size_invitro))
a+geom_smooth(method = 'lm', alpha = 0.3, colour = "blue4", fill='lightskyblue',linetype=1)+
  geom_point(size = 2.5, colour = "grey25", alpha = 0.8)+
  stat_cor(method = "pearson")+
  theme_classic()+
  labs(title = "Effect Size conc lipids Exp#3 vs. GPMV" , 
       x = "GPMV", 
       y = "Primary culture #3")+
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
ggsave(filename =  "./Figure_Panels/Fig.3b.pdf", height = 5, width = 5, useDingbats=FALSE)
