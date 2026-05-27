
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)

source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/EffectSize_DESI_annotated_lipids.Rdata")
load("./Output_Data/Lpd.MsD.Invivo_Age_ES.w.RT.Rdata")
load("./Output_Data/Lpd.MsD.Invitro_Qui_Age_ES.w.RT.Rdata")

DESI <- ES.lpd.anno.inDESI

Invitro.es <- MsD.lpd.rmv.abc(Invitro.Q.lpd.ES.w.RT) %>% 
  ungroup() %>% 
  rowwise() %>% 
  filter(LipidIon %in% DESI$LipidIon) %>% 
  relocate(es_g, .after = "LipidIon") %>% 
  rename("es_g.Exp3" = "es_g") %>% 
  relocate(es_g.Exp3, .after = "LipidIon") %>% 
  filter(!(LipidIon == "PE(P-18:1_18:1)" & Average.Rt.min. == 13.351)) %>% 
  filter(!(LipidIon == "PI(18:0_20:3)" & Average.Rt.min. == 11.519))

Invivo.es <- MsD.lpd.rmv.abc(Invivo.lpd.ES.w.RT) %>%
  ungroup() %>% 
  rowwise() %>% 
  filter(LipidIon %in% DESI$LipidIon) %>%  
  rename("es_g.invivo" = "es_g")

all.lc <- list(Invitro.es, Invivo.es) %>% 
  reduce(full_join, by = "LipidIon") #11

all.lc.DESI.org <- left_join(all.lc, DESI, by = "LipidIon") %>%  
  ungroup() %>% 
  select(LipidIon, starts_with("es")) %>% 
  pivot_longer(-LipidIon, names_to = "Exp", values_to = "Exp_es") %>% 
  mutate(Experiment = case_when(
    grepl("Exp3", Exp) ~ "In vitro",
    grepl("invivo", Exp) ~ "In vivo",
    Exp == "es" ~ "DESI"
  ))

ES.mean <- all.lc.DESI.org %>% 
  group_by(LipidIon) %>% 
  summarise(MeanES = mean(Exp_es, na.rm = TRUE), SEM = sd(Exp_es, na.rm = TRUE)/sqrt(sum(!is.na(Exp_es)))) 


CI95 <- function(df, mean, sem){
  df <- df %>% 
    rowwise() %>% 
    mutate(`Margin of Error` = 1.96 * {{sem}}) %>% 
    mutate(CI.lower = {{mean}} -`Margin of Error`) %>% 
    mutate(CI.upper = {{mean}} +`Margin of Error`)
}

ES.mean.CI.LoUp <- CI95(ES.mean, MeanES, SEM)

all.ES.p <- left_join(all.lc.DESI.org, ES.mean.CI.LoUp, by = "LipidIon") %>% 
  mutate(Sig = ifelse( CI.lower > 0 & CI.upper > 0 | CI.lower < 0 & CI.upper < 0,
                       "Significant", "Not significant"))

all.ES.p.org <- all.ES.p %>% 
  mutate(Cat = ifelse(Experiment == "DESI", "DESI", "LC-MS"))

all.ES.p.org$Experiment <- factor(all.ES.p.org$Experiment, levels = c("In vitro", "In vivo", "DESI"))

sig.ES <- all.ES.p.org %>% 
  filter(Sig == "Significant") #6

nosig.ES <- all.ES.p.org %>% 
  filter(!Sig == "Significant") #27

lpd.DESI.all <- bind_rows(sig.ES, nosig.ES) %>% 
  mutate(LipidIon.l = ifelse(Sig == "Significant", 
                             paste0(LipidIon, "<b>*</b>"), LipidIon))


a <-ggplot(lpd.DESI.all, aes(x = fct_reorder(LipidIon.l, MeanES), y = Exp_es))
a+
  # geom_point(aes(shape = Experiment, color = Cat), alpha = 0.8, size = 3) +
  geom_point(aes(shape = Experiment, color = Cat), alpha = 0.85, size = 3.5) +
  scale_shape_manual(values = c(17, 25, 18)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM), alpha = 0.75, width = 0.2, colour = "grey15",) +
  stat_summary(aes(x=LipidIon.l,y=MeanES), fun=mean, geom = "point", size=6, shape=20, alpha = 0.75, colour = "grey15",) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("seagreen", "tan2"))+
  coord_flip() + 
  labs(title = "DESI vs. In vitro and In vivo LC-MS", x = "", y = "Effect size (Old vs. Young)") +
  theme(axis.text.y = element_markdown(colour = "black", size = 6))
ggsave(filename = "./Figure_Panels/Fig.3D.pdf", width = 5, height = 5, useDingbats = FALSE)

## ==== highlight Mboat2 responsive lipids====
load("./Output_Data/Mboat2.responsive_lipid_list.Rdata")

DESI.ovlp.ls <- unique(lpd.DESI.all$LipidIon)[unique(lpd.DESI.all$LipidIon) %in% res.lpd.ls]
DESI.ovlp.ls
# character(0)
