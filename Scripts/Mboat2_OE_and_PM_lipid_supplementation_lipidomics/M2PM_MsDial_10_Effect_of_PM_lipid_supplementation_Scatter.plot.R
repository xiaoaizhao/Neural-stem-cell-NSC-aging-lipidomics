##Prerequisite: Need to run the following 3 scripts first
### 'Fig1.Lipid_2.invitro.studies.R'
### 'Fig.4.OE_Lipids_responsive_to_Mboat2_manipulation.R'
### 'Fig4.KO_Lipids_responsive_to_Mboat2_manipulation.R'
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

## ==== highlight Mboat2 responsive lipids====
load("./Output_Data/OE_responsive.lpd.df.Rdata")
load("./Output_Data/KO_responsive.lpd.df.Rdata")

res.lpd.ls.all <- bind_rows(padj.sig.KO.lpd., padj.sig.OE.lpd.) %>% 
  select(LipidIon, Aging.dir)

PM.supp.effect.lb.p <- left_join(PM.supp.effect.lb, res.lpd.ls.all, by = "LipidIon") %>% 
  mutate(Label.all = ifelse(is.na(Aging.dir), "", LipidIon)) %>% 
  mutate(Aging.D = ifelse(is.na(Aging.dir), "Not sig", Aging.dir))

PM.supp.effect.lb.p$Aging.D <- factor(PM.supp.effect.lb.p$Aging.D, levels = c("Not sig", "Higher in young", "Higher in old"))

a <- ggplot(PM.supp.effect.lb.p, aes(`YGPMV treatment in cell`, `OvY GPMV lipids`))
a+
  geom_point(aes(color = Aging.D), size = 1, alpha = 0.75)+
  scale_color_manual(values = c("grey39", "cyan3", "magenta3"))+
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
                  max.overlaps = 45)
ggsave(filename = paste0("./Figure_Panels/Fig.5b.pdf"), width = 6, height = 4, useDingbats = FALSE)

