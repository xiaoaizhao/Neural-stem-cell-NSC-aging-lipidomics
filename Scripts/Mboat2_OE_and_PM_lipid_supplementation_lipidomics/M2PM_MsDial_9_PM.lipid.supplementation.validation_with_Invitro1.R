### Validation on young PM lipid supplementation by lipidomics

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggthemes)
library(ggpubr)

source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")
##==== Effect size between PM lipid supplemented cells and control cells ====
load(file = "./Output_Data/PM_lipid_supplementation_by_treatment_all_ages.Rdata")

##==== GPMV data ====
load(file = "./Output_Data/GPMV_MsD.Raw.conc.Norm_Impt_291_lipids.Rdata")
GPMV <- GPMV.impt.norm.raw.conc.lpd.MsD %>% 
  rownames_to_column(., var = "LipidIon")
  
GPMV.df <- MsD.lpd.rmv.abc(GPMV)

##==== In vitro data ====
load("./Output_Data/Invitro_conc.lipid.qNSC.Rdata")

Invitro <- MsD.lpd.rmv.abc(Invitro.q.conc)

GPMV.vs.Cell <- inner_join(GPMV.df, Invitro, by = "LipidIon") %>% 
  pivot_longer(-LipidIon, names_to = "Samples", values_to = "Concentration") %>% 
  mutate(Condition = ifelse(grepl("_qNSC-Q", Samples), "Control", "Treatment"))

##==== Effect size between whole cell and GPMV dataset ====
es.GPMVvCell <- es.g.treat.func(GPMV.vs.Cell, LipidIon, Condition, Concentration, Samples)  # positive means higher in GPMV

es.PM.supp <- MsD.lpd.rmv.abc(PM.supp.ef.sz)

es.supp.val <- inner_join(es.PM.supp, es.GPMVvCell, by = c("LipidIon")) %>% 
  rename( "Effect_size_GPMV.sup" = "es_g.x") %>% 
  rename( "Effect_size_GPMV.vs.wholecell" = "es_g.y")

es.supp.val$Age <- factor(es.supp.val$Age, levels = c("Young", "Old"))

sup.val.p <- es.supp.val %>% 
  mutate(Sup.cat = ifelse(Effect_size_GPMV.sup > 0, "Higher w/supp", "Lower w/supp")) %>% 
  mutate(Overall.cat = ifelse(Sup.cat == "Higher w/supp",
         case_when(
           Effect_size_GPMV.vs.wholecell > 0 ~ "Higher w/GPMV",
           Effect_size_GPMV.vs.wholecell < 0 ~ "Lower w/GPMV",
         ), "Lower w/supp"))

sup.val.p$Sup.cat <- factor(sup.val.p$Sup.cat, levels = c( "Lower w/supp", "Higher w/supp"))
sup.val.p$Overall.cat <- factor(sup.val.p$Overall.cat, levels = c( "Lower w/supp", "Higher w/GPMV", "Lower w/GPMV"))
sup.val.p <- sup.val.p %>% 
  arrange(Overall.cat)

pal2 <- c( "grey80", "#737070", "#008585")

a <- ggplot(sup.val.p, aes(reorder(LipidIon, Effect_size_GPMV.vs.wholecell), Effect_size_GPMV.vs.wholecell))
a+ geom_point(aes(color = Overall.cat), alpha = 0.95, size = 3)+
  theme_classic()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(title = "Lipids higher in cells with supplementation")+
  scale_x_discrete(expand = c(0.1,0.1)) +
  scale_color_manual(values = pal2) +
  theme(legend.position = "bottom") +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Effect size GPMV vs. Whole cell") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~Age)
ggsave(filename = "./Figure_Panels/fig.S14C.pdf", width = 6, height = 4)

#### ==== add pie chart to indicate what percentage of lipids are more abundant in GPMV vs. whole cell====####
pie.df <- sup.val.p %>% 
  mutate(GPMV.cat = ifelse(Effect_size_GPMV.vs.wholecell > 0, "Higher w/GPMV", "Lower w/GPMV")) %>% 
  filter(Sup.cat == "Higher w/supp") %>% 
  group_by(Age, GPMV.cat, Sup.cat) %>% 
  tally()

sup.lipid.sum <- sup.val.p %>% 
  mutate(GPMV.cat = ifelse(Effect_size_GPMV.vs.wholecell > 0, "Higher w/GPMV", "Lower w/GPMV")) %>% 
  filter(Sup.cat == "Higher w/supp") %>% 
  group_by(Age) %>% 
  summarise(nLipid = n())

pie.df.p <- left_join(pie.df, sup.lipid.sum, by = "Age") %>% 
  mutate(pctLipid = (n/nLipid)*100)

pie.df.p$GPMV.cat <- factor(pie.df.p$GPMV.cat, levels = c("Higher w/GPMV", "Lower w/GPMV"))
pie.df.p$Age <- factor(pie.df.p$Age, levels = c("Young", "Old"))
pie.df.p <- pie.df.p %>% 
  arrange(Age, GPMV.cat) %>% 
  mutate(ypos = cumsum(pctLipid) - 0.5 * pctLipid)


pal.p <- c(  "#737070", "#008585")
a <- ggplot(pie.df.p, aes(x="", y = pctLipid, fill = GPMV.cat), alpha = 1)
a+  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = ypos, label = paste0(round(pctLipid,1), "%")),
            color = "white", size = 4) +
  theme_classic()+
  scale_fill_manual(values = pal.p) +
  facet_wrap(~Age)+
  theme_void()+
  theme(legend.position = "right")
ggsave(filename = "./Figure_Panels/fig.S14C.pie.pdf",
       width = 8, height = 4, useDingbats = FALSE)
