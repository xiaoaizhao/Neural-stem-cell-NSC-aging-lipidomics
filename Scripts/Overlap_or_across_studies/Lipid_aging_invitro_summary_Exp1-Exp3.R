
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)

source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_Data/Ef_Size_Lipid_Age_qNSC_InVitro_LC.Rdata") # Primary culture #1
load("./Output_Data/Ef_Size_CONC.Lipid_Exp2_all_KO.Rdata") # Primary culture #2
load("./Output_Data/Exp3_Qui_Age_ES.conc.lipids.Rdata")# Primary culture #3

Exp2.ctrl.conc <- Exp2.CONC.lpd.es.g.allKO %>% 
  filter(KO == "N") 

CI95 <- function(df, mean, sem){
  df <- df %>% 
    rowwise() %>% 
    mutate(`Margin of Error` = 1.96 * {{sem}}) %>% 
    mutate(CI.lower = {{mean}} -`Margin of Error`) %>% 
    mutate(CI.upper = {{mean}} +`Margin of Error`)
}

E3 <- E3.Q.conc.AG.ES
E2 <- ether.rename(Exp2.ctrl.conc)
E1 <- ether.rename(InVitro.lpd.es.g)

c2.in.3lc<- list(E1, 
                 E2,
                 E3) %>% 
  reduce(full_join, by = c("LipidIon")) #949
  
c2.lpd <- c2.in.3lc %>% 
  select(LipidIon, starts_with("es_g")) %>% 
  mutate(Detected = 3 - rowSums(is.na(across(starts_with("es_g"))))) %>% 
  filter(Detected >=2) #301 lipids


c2.in.3lc.lpd <- c2.in.3lc %>% 
  filter(LipidIon %in% c2.lpd$LipidIon) %>% 
  select(LipidIon, starts_with("es_g")) %>% 
  pivot_longer(-LipidIon, names_to = "Exp", values_to = "Effect size") %>% 
  group_by(LipidIon) %>% 
  summarise(MeanES = mean(`Effect size`, na.rm = TRUE), SEM = sd(`Effect size`, na.rm = TRUE)/sqrt(sum(!is.na(`Effect size`)))) 

c2.in.3lc.lpd.CI.LoUp <- CI95(c2.in.3lc.lpd, MeanES, SEM)

E1.p <-  E1 %>% 
  filter(LipidIon %in% c2.in.3lc.lpd.CI.LoUp$LipidIon) %>% 
  select(LipidIon, es_g) %>% 
  mutate(Exp = "Primary culture #1")

E2.p <-E2 %>% 
  filter(LipidIon %in% c2.in.3lc.lpd.CI.LoUp$LipidIon) %>% 
  select(LipidIon, es_g) %>% 
  mutate(Exp = "Primary culture #2")

E3.p <- E3 %>% 
  filter(LipidIon %in% c2.in.3lc.lpd.CI.LoUp$LipidIon) %>% 
  select(LipidIon, es_g) %>% 
  mutate(Exp = "Primary culture #3")


c2.in.3lc.df <- bind_rows(E1.p, E2.p, E3.p)
p.c2.in.3lc.df <- left_join(c2.in.3lc.df, c2.in.3lc.lpd.CI.LoUp, by = "LipidIon") %>% 
   mutate(Sig = ifelse( CI.lower > 0 & CI.upper > 0 | CI.lower < 0 & CI.upper < 0,
         "Significant", "Not significant")) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
    .x %>% 
      mutate(star.pos = case_when(
        Sig == "Significant" & MeanES > 0 ~  max(es_g) + 0.8, 
        Sig == "Significant" & MeanES < 0 ~ min(es_g) - 0.8)
) 
  }) 
 
p.c2.in.3lc.df$Exp <- factor(p.c2.in.3lc.df$Exp, levels = c( "Primary culture #1", "Primary culture #2", "Primary culture #3"))


save(p.c2.in.3lc.df, file = "./Output_Data/Lipids.2outof3LC.invitro.allfeatures.Rdata")


lpd.c2.in.3lc.df.sig <- p.c2.in.3lc.df %>% 
  filter(Sig == "Significant")

save(lpd.c2.in.3lc.df.sig, file = "./Output_Data/Lipids.2outof3LC.invitro.Sig.features.Rdata")

load("./Output_Data/Lipids.2outof3LC.invitro.Sig.features.Rdata")
a <-ggplot(lpd.c2.in.3lc.df.sig, aes(x = fct_reorder(LipidIon, MeanES), y = es_g))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 2, fill = "gold") +
  scale_shape_manual(values = c(24, 2, 17, 10)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM), colour = "grey15", alpha = 0.75, width = 0.2) +
  stat_summary(aes(x=LipidIon,y=MeanES), fun=mean, geom = "point", size=2, shape=20, alpha = 0.75, colour = "grey15") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 4), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() + 
  labs(title = "Significant 2 out of 3 LC in vitro lipids", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = lpd.c2.in.3lc.df.sig %>% 
           filter(Sig == "Significant"), 
           aes(x = fct_reorder(LipidIon, MeanES), y = star.pos),
         pch=8, 
         size=0.8, stroke = 0.7, alpha = 0.75,
         colour="black") +
  theme(legend.position = "bottom")
ggsave(filename = "./Figure_Panels/EDFig.1h.pdf", width = 4, height = 12, useDingbats = FALSE)

## ====Top 30% biggest effect size to plot in main figure ====
lpd.c2.in.3lc.df.sig.top30.df <- lpd.c2.in.3lc.df.sig %>% 
  group_by(LipidIon) %>% 
  summarise(MeanES = unique(MeanES)) %>% 
  mutate(percentile = percent_rank(abs(MeanES))) %>% 
  filter(percentile > 0.7) #38 lipids

lpd.c2.in.3lc.df.sig.top30.p <- lpd.c2.in.3lc.df.sig %>% 
  filter(LipidIon %in% lpd.c2.in.3lc.df.sig.top30.df$LipidIon)

a <-ggplot(lpd.c2.in.3lc.df.sig.top30.p, aes(x = fct_reorder(LipidIon, MeanES), y = es_g))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 2, fill = "gold",) +
  scale_shape_manual(values = c(24, 2, 17, 10)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM), colour = "grey15", alpha = 0.75, width = 0.2) +
  stat_summary(aes(x=LipidIon,y=MeanES), fun=mean, geom = "point", size=2, shape=20, alpha = 0.75, colour = "grey15") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 4), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() + 
  labs(title = "Significant 2 out of 3 LC in vitro lipids - top 30%", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = lpd.c2.in.3lc.df.sig.top30.p %>% 
           filter(Sig == "Significant"), 
           aes(x = fct_reorder(LipidIon, MeanES), y = star.pos),
         pch=8, 
         size=0.8, stroke = 0.7, alpha = 0.75,
         colour="black") +
  theme(legend.position = "bottom")
ggsave(filename = "./Figure_Panels/Fig.1f.pdf", width = 4, height = 5, useDingbats = FALSE)
