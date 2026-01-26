### Plot effect size by intervention on highlighted lipids

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(patchwork)
source("./Scripts/Function_scripts/Pre-processing_functions.R")
## ==== Effect on Mboat2 KO ====
### ==== old ====
# In this dataset, positive effect size means higher in control, negative effect size means higher in Mboat2 KO
load("./Output_Data/Exp2_lpd.es.by.KO.in.old.Rdata")

up.ls <- c("LPC(O-16:0)", "LPC(O-18:1)","PI(18:0_20:4)", "PE(P-18:1_18:1)")
down.ls <- c("PI(18:0_20:3)", "PC(16:0_20:3)")

M2KO.old.up <- MsD.lpd.rmv.abc(all.ko.es.old) %>% 
  filter(KO == "M" & LipidIon %in% up.ls) %>% 
  mutate(Exp = "Mboat2 KO")

M2KO.old.dw <- MsD.lpd.rmv.abc(all.ko.es.old) %>% 
  filter(KO == "M" & LipidIon %in% down.ls) %>% 
  mutate(Exp = "Mboat2 KO")


## ==== Effect on Mboat2 OE ====
### ==== old ====
# In this dataset, positive effect size means higher in Mboat2 OE, negative effect size means higher in control
load("./Output_Data/M2OE_lpd.es.by.OE.in.Old.Rdata")

up.ls <- c("LPC(O-16:0)", "LPC(O-18:1)","PI(18:0_20:4)", "PE(P-18:1_18:1)")
down.ls <- c("PI(18:0_20:3)", "PC(16:0_20:3)")

M2OE.old.up <- M2.old.OE.ES%>% 
  filter( LipidIon %in% up.ls) %>% 
  mutate(Exp = "Mboat2 OE")

M2OE.old.dw <- M2.old.OE.ES%>% 
  filter( LipidIon %in% down.ls) %>% 
  mutate(Exp = "Mboat2 OE")

## ==== all interventions together====
### old samples only
#### upwith age
up.age.all <- bind_rows(M2KO.old.up, M2OE.old.up)
up.age.all$Exp <- factor(up.age.all$Exp, levels = c("Mboat2 OE", "Mboat2 KO"))


#### down with age
down.age.all <- bind_rows(M2KO.old.dw, M2OE.old.dw)
down.age.all$Exp <- factor(down.age.all$Exp, levels = c( "Mboat2 OE", "Mboat2 KO"))


#### ==== make single plot for LPC(O-16:0) and PI(18:0_20:3)
LPC.up <- up.age.all %>% 
  filter(LipidIon == "LPC(O-16:0)")

up.p <-ggplot(LPC.up, aes(x = Exp, y = es_g))+
  geom_point(aes(shape = Exp), colour = "#CC66CC", alpha = 0.95, size = 3) +
  scale_shape_manual(values = c(7,11, 14)) + 
  geom_errorbar(aes(ymin = es_g - se_g, ymax = es_g + se_g), colour = "black", alpha = 0.75, width = 0.05) +
  # stat_summary(aes(x=LipidIon,y=MeanES), fun=mean, geom = "point", size=2, shape=20, alpha = 0.75, colour = "grey15") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 7.5), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "LPC(O-16:0) increases with age", x = "", y = "Effect size (Treatment vs. Control)") +
  theme(legend.position = "none") +
  coord_flip()

PI.dw <- down.age.all %>% 
  filter(LipidIon == "PI(18:0_20:3)") %>% 
  filter(!(es_g > 0 & Exp == "Mboat2 KO")) %>% 
  filter(!(es_g < 0 & Exp == "Young GPMV supp"))  

down.p <-ggplot(PI.dw, aes(x = Exp, y = es_g))+
  geom_point(aes(shape = Exp), colour = "#00CDCD", alpha = 0.95, size = 3) +
  scale_shape_manual(values = c( 7,11, 14)) + 
  geom_errorbar(aes(ymin = es_g - se_g, ymax = es_g + se_g), colour = "black", alpha = 0.75, width = 0.05) +
  # stat_summary(aes(x=LipidIon,y=MeanES), fun=mean, geom = "point", size=2, shape=20, alpha = 0.75, colour = "grey15") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 7.5), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "PI(18:0_20:3) decreases with age", x = "", y = "Effect size (Treatment vs. Control)") +
  theme(legend.position = "none") +
  coord_flip()

up.p / down.p
ggsave(filename = "./Figure_Panels/EDFig.11f.pdf", height = 4, width = 2, useDingbats = FALSE)
