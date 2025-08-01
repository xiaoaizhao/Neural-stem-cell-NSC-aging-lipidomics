## Validation on young PM lipid supplementation by lipidomics

##Prerequisite: Need to run scripts #0 - #5 in GPMV lipidomics folder and scripts #0 - #5 in Primary culture #2 lipidomics first

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggthemes)
library(ggpubr)

source("./Scripts/Function_scripts/Effect_size_functions.R")
##==== Effect size between PM lipid supplemented cells and control cells ====
load(file = "./Output_Data/PM_lipid_supplementation_by_treatment_all_ages.Rdata")

##==== GPMV data ====
load(file = "./Output_Data/GPMV_Norm_Impt_log2_conc_483_lipids.Rdata")
conc.gpmv <- 2^GPMV_Impt_norm_conc_all %>% 
  rownames_to_column(., var = "LipidIon")
  
GPMV.df <- ether.rename(conc.gpmv)

##==== Primary culture #2 ====
load("./Output_Data/Exp2_Norm_Impt_backtoraw_conc612_lipids.Rdata")

Exp2.ctrl <- raw_conc.exp2 %>% 
  select(contains("_N")) %>% 
  rownames_to_column(var = "LipidIon")

Exp2.ctrl.df <- ether.rename(Exp2.ctrl)

GPMV.vs.Exp2Cell <- inner_join(GPMV.df, Exp2.ctrl.df, by = "LipidIon") %>% 
  pivot_longer(-LipidIon, names_to = "Samples", values_to = "Concentration") %>% 
  mutate(Condition = ifelse(grepl("_N", Samples), "Control", "Treatment"))

##==== Effect size between whole cell and GPMV dataset ====
es.GPMVvCell <- es.g.treat.func(GPMV.vs.Exp2Cell, LipidIon, Condition, Concentration, Samples)  # positive means higher in GPMV

es.PM.supp <- ether.rename(PM.supp.ef.sz)

es.supp.val <- inner_join(es.PM.supp, es.GPMVvCell, by = c("LipidIon")) %>% 
  rename( "Effect_size_GPMV.sup" = "es_g.x") %>% 
  rename( "Effect_size_GPMV.vs.wholecell" = "es_g.y")

es.supp.val$Age <- factor(es.supp.val$Age, levels = c("Young", "Old"))

sup.val.p <- es.supp.val %>% 
  mutate(Sup.cat = ifelse(Effect_size_GPMV.sup > 0, "Higher w/supp", "Lower w/supp"))

sup.val.p$Sup.cat <- factor(sup.val.p$Sup.cat, levels = c( "Lower w/supp", "Higher w/supp"))

sup.val.p <- sup.val.p %>% 
  arrange(Sup.cat)

pal2 <- c( "grey80", "#330000")

a <- ggplot(sup.val.p, aes(reorder(LipidIon, Effect_size_GPMV.vs.wholecell), Effect_size_GPMV.vs.wholecell))
a+ geom_point(aes(color = Sup.cat), alpha = 0.75, size = 3)+
  theme_classic()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(title = "Lipids higher in cells with supplementation")+
  scale_x_discrete(expand = c(0.1,0.1)) +
  scale_color_manual(values = pal2) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Effect size GPMV vs. Whole cell") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~Age)
ggsave(filename = "./Figure_Panels/EDFig.13c.pdf", width = 6, height = 4)

