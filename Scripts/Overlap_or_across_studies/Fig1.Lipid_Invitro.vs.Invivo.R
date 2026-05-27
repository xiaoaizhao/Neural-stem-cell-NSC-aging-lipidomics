

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(ggtext)
source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Lpd.MsD.Invivo_Age_ES.w.RT.Rdata")
load("./Output_Data/Lpd.MsD.Invitro_Qui_Age_ES.w.RT.Rdata")

CI95 <- function(df, mean, sem){
  df <- df %>% 
    rowwise() %>% 
    mutate(`Margin of Error` = 1.96 * {{sem}}) %>% 
    mutate(CI.lower = {{mean}} -`Margin of Error`) %>% 
    mutate(CI.upper = {{mean}} +`Margin of Error`)
}

Invivo <- MsD.lpd.rmv.abc(Invivo.lpd.ES.w.RT) %>% 
  select(LipidIon, Average.Rt.min., es_g) %>% 
  rename("es_g.invivo" = "es_g")

Invitro <- MsD.lpd.rmv.abc(Invitro.Q.lpd.ES.w.RT) %>% 
  select(LipidIon, Average.Rt.min., es_g) %>% 
  rename("es_g.Exp3" = "es_g")

Invitro.Invivo.lpd <- inner_join(Invitro, Invivo, by = "LipidIon") %>% 
  mutate(DeltaRT = Average.Rt.min..y - Average.Rt.min..x)
  

lpd.Invitro.Invivo.l <- Invitro.Invivo.lpd %>% 
  filter(!(LipidIon == "LPE(18:0)" & Average.Rt.min..x == 3.990)) %>% 
  ungroup() %>% 
  select(LipidIon, starts_with("es")) %>% 
  pivot_longer(-LipidIon, names_to = "Exp", values_to = "Effect size")

sum.Invitro.Invivo <- lpd.Invitro.Invivo.l %>% 
  group_by(LipidIon) %>% 
  summarise(MeanES = mean(`Effect size`, na.rm = TRUE), SEM = sd(`Effect size`, na.rm = TRUE)/sqrt(sum(!is.na(`Effect size`)))) 

lpd.Invitro.Invivo.CI.LoUp <- CI95(sum.Invitro.Invivo, MeanES, SEM)

c2of2sig.p <- inner_join(lpd.Invitro.Invivo.l, lpd.Invitro.Invivo.CI.LoUp, by = "LipidIon") %>% 
  mutate(Sig = ifelse( CI.lower > 0 & CI.upper > 0 | CI.lower < 0 & CI.upper < 0,
                       "Significant", "Not significant")) %>% 
  group_by(LipidIon) %>% 
  mutate(Exp = case_when(
    grepl("Exp3", Exp) ~ "In vitro",
    grepl("invivo", Exp) ~ "In vivo",
  )) %>% 
  rename("es_g" = "Effect size") %>% 
  mutate(LipidIon.l = ifelse(Sig == "Significant", 
                             paste0(LipidIon, "<b>*</b>"), LipidIon))

c2of2sig.p$Exp <- factor(c2of2sig.p$Exp, levels = c( "In vitro", "In vivo"))

sig.invivo <- c2of2sig.p %>% 
  filter(Sig == "Significant") #2 only 2 significant lipids
  
## ==== Calculate p value based on CI ====


a <-ggplot(c2of2sig.p, aes(x = fct_reorder(LipidIon.l, MeanES), y = es_g))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 3.5) +
  scale_shape_manual(values = c(17, 25)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM), colour = "grey15", alpha = 0.75, width = 0.2) +
  stat_summary(aes(x=LipidIon.l,y=MeanES), fun=mean, geom = "point", size=3.5, shape=20, alpha = 0.75, colour = "grey15") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 8), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # scale_color_manual(values = c("cornflowerblue", "yellowgreen", "firebrick3"))+
  coord_flip() + 
  labs(title = "In vitro vs. In vivo", x = "", y = "Effect size (Old vs. Young)") +
  theme(axis.text.y = element_markdown(colour = "black", size = 6)) +
  theme(legend.position = "bottom")
ggsave(filename = "./Figure_Panels/Fig.2B.pdf", width = 5, height = 5, useDingbats = FALSE)

## ==== highlight Mboat2 responsive lipids====
load("./Output_Data/Mboat2.responsive_lipid_list.Rdata")

invivo.ovlp.ls <- unique(c2of2sig.p$LipidIon)[unique(c2of2sig.p$LipidIon) %in% res.lpd.ls]
invivo.ovlp.ls
