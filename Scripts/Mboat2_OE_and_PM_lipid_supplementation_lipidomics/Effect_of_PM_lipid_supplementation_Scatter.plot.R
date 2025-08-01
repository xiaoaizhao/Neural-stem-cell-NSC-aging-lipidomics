##Prerequisite: Need to run `Effect_of_PM_lipid_supplementation_Pie.chart.R` first
library(tidyverse)
library(ggrepel)
library(ggpubr)
setwd(rstudioapi::getActiveProject())
rm(list = ls())

# upload dataframe, top hit column T/F indicates whether lipids show the top 30 biggest effect size change on both axis
load("./Output_Data/PM.lipid.supp.effect.by.GPMV.aging.ES.Rdata")

## ====isolate examples lipids ====####
PM.supp.effect.lb <- PM.supp.effect %>%
  mutate(YO.gpmv.rk.raw = percent_rank(`OvY GPMV lipids`)) %>%
  mutate(gpmv.supp.rk.raw = percent_rank(`YGPMV treatment in cell`)) %>%
  group_by(Age) %>%
  group_modify(~{
    .x %>%
      mutate(old.rmv = ifelse(
        YO.gpmv.rk.raw>0.7 & gpmv.supp.rk.raw<0.3, "T", "F"
      )) %>%
      mutate(yng.accu = ifelse(
        YO.gpmv.rk.raw<0.3 & gpmv.supp.rk.raw>0.7, "T", "F"
      ))
  }) %>%
  mutate(old.rmv.lb = ifelse(old.rmv == "T", LipidIon, "")) %>%
  mutate(yng.accu.lb = ifelse(yng.accu == "T", LipidIon, ""))


## ==== In one plot - increased lipids with supplementation (associated with young GPMV) AND decreased lipids with supplementation (associated with old GPMV)====
PM.supp.effect.lb.p <- PM.supp.effect.lb %>% 
  mutate(Label.all = paste0(yng.accu.lb, old.rmv.lb))

a <- ggplot(PM.supp.effect.lb.p, aes(`YGPMV treatment in cell`, `OvY GPMV lipids`))
a+
  geom_point(colour = "grey39", size = 1, alpha = 0.75)+
  theme_classic()+
  labs(title = "Effect of PM lipid supplementation" , 
       y = "Effect size - (old vs. young GPMV)", 
       x = "Effect size - (young GPMV supplementation vs. control)")+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(legend.position = "none") + 
  facet_wrap(~Age, scales = "fixed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(aes(label = Label.all), fontface = 'plain',
                  size = 1.5,colour = "black",
                  box.padding = unit(0.55, "lines"),
                  seed = 1234,
                  min.segment.length = 0,
                  max.overlaps = 25)
ggsave(filename = "./Figure_Panels/Fig.5b.pdf", width = 6, height = 4, useDingbats = FALSE)
