### Visualize the effect of young PM lipid supplementation by category in pie chart
### Based on the effect size comparison between young and old GPMV, as well as effect size between cells with or without young PM lipid supplementation


setwd(rstudioapi::getActiveProject())
rm(list = ls())

library(tidyverse)
library(ggpubr)
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load(file = "./Output_Data/PM_lipid_supplementation_by_treatment_all_ages.Rdata")

# Effect size of old GPMV lipids vs. young GPMV lipids (lipids with quantitative standards only)
load("./Output_Data/Lpd.MsD.GPMV_Age_ES.Rdata")
ovy.gpmv <- MsD.lpd.rmv.abc(GPMV.lpd.ES.w.RT) %>% 
  select(LipidIon, es_g)

# Effect size of old cells w/Y GPMV vs. old cells w/Control
wYGPMV <- MsD.lpd.rmv.abc(PM.supp.ef.sz) %>% 
  select(LipidIon, es_g, Age)

# Merge to data frame and correlation analysis
PM.supp.effect <- inner_join(ovy.gpmv, wYGPMV, by = "LipidIon") %>% 
  ungroup() %>% 
  rename("OvY GPMV lipids" = "es_g.x") %>% 
  rename("YGPMV treatment in cell" = "es_g.y") %>% 
  mutate(Cat.GPMV = ifelse(`OvY GPMV lipids` < 0, "Enriched in Y GPMV", "Enriched in O GPMV")) %>% 
  mutate(Cat.GPMV.supp = ifelse(`YGPMV treatment in cell` < 0, "Depleted with Y GPMV supp", "Enriched with Y GPMV supp")) %>% 
  mutate(YO.gpmv.rk = percent_rank(abs(`OvY GPMV lipids`))) %>% 
  mutate(gpmv.supp.rk = percent_rank(abs(`YGPMV treatment in cell`))) %>% 
  mutate(Top.hit = ifelse(
    YO.gpmv.rk > 0.3 | gpmv.supp.rk > 0.3,
    "T", "F"
  ))

PM.supp.effect$Cat.GPMV <- factor(PM.supp.effect$Cat.GPMV, levels = c("Enriched in Y GPMV", "Enriched in O GPMV"))
PM.supp.effect$Cat.GPMV.supp <- factor(PM.supp.effect$Cat.GPMV.supp, levels = c("Depleted with Y GPMV supp", "Enriched with Y GPMV supp"))

PM.supp.effect$Age <- factor(PM.supp.effect$Age, levels = c("Young", "Old"))
save(PM.supp.effect, file = "./Output_Data/PM.lipid.supp.effect.by.GPMV.aging.ES.Rdata")

#### ==== pie chart to tally lipids from each quadrant====####
tally.y.o <- PM.supp.effect %>% 
  group_by(Age, Cat.GPMV.supp, Cat.GPMV) %>% 
  summarise(nLipid = n()) %>% 
  mutate(grp = paste0(Cat.GPMV.supp, "_", Cat.GPMV)) %>% 
  ungroup(Cat.GPMV.supp, Cat.GPMV) %>% 
  group_modify(~{
    .x %>% 
      mutate(pctLipid = nLipid/sum(nLipid)*100)
  }) 

tally.y.o$Cat.GPMV.supp <- factor(tally.y.o$Cat.GPMV.supp, levels = c("Enriched with Y GPMV supp",
                                                                      "Depleted with Y GPMV supp"))

tally.y.o <- tally.y.o %>% 
  arrange(Age, Cat.GPMV.supp) %>% 
  mutate(ypos = cumsum(pctLipid) - 0.5 * pctLipid)
tally.y.o$Age <- factor(tally.y.o$Age, levels = c("Young", "Old"))

pal.p <- c( "#0e18f0", "#0e89f0", "#f0760e", "#f0af0e")
a <- ggplot(tally.y.o, aes(x="", y = pctLipid, fill = grp), alpha = 0.7)
a+  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = ypos, label = paste0(round(pctLipid,1), "%")),
            color = "white", size = 4) +
  theme_classic()+
  scale_fill_manual(values = pal.p) +
  facet_wrap(~Age)+
  theme_void()+
  theme(legend.position = "right")
ggsave(filename = "./Figure_Panels/EDFig.13f.pdf",
       width = 8, height = 4, useDingbats = FALSE)
