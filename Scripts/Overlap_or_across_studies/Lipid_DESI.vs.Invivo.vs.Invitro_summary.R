
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)

source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_Data/EffectSize_DESI_annotated_lipids.Rdata")
load("./Output_Data/Ef_Size_CONC.Lipid_InVivo.Rdata")
load("./Output_Data/Lipids.2outof3LC.invitro.allfeatures.Rdata")

c2.3.all <- p.c2.in.3lc.df %>% 
  select(LipidIon, MeanES) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
    .x %>% 
      summarise(SumES = unique(MeanES))
  })

invitro <- c2.3.all %>% 
  filter(LipidIon %in% ES.lpd.anno.inDESI$LipidIon) #12
  
Invivo.es <- ether.rename(Invivo.CONC.lpd.es.g) %>% 
  filter(LipidIon %in% ES.lpd.anno.inDESI$LipidIon)

all.lc <- list(invitro, Invivo.es) %>% 
  reduce(full_join, by = "LipidIon")

all.lc.DESI.org <- left_join(all.lc, ES.lpd.anno.inDESI, by = "LipidIon") %>%  
  select(LipidIon, SumES, es_g, es) %>% 
  pivot_longer(-LipidIon, names_to = "Exp", values_to = "Exp_es") %>% 
  mutate(Experiment = case_when(
    Exp == "SumES" ~ "Invitro summary",
    Exp == "es_g" ~ "In vivo",
    Exp == "es" ~ "DESI",
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
         "Significant", "Not significant")) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
    .x %>% 
      mutate(star.pos = case_when(
        Sig == "Significant" & MeanES > 0 ~  max(Exp_es, na.rm = TRUE) + 0.8, 
        Sig == "Significant" & MeanES < 0 ~ min(Exp_es, na.rm = TRUE) - 0.8)
)
  })

all.ES.p.org <- all.ES.p %>% 
  mutate(Cat = ifelse(Experiment == "DESI", "DESI", "LC-MS")) 

all.ES.p.org$Experiment <- factor(all.ES.p.org$Experiment, levels = c("Invitro summary", "In vivo", "DESI"))


a <-ggplot(all.ES.p.org, aes(x = fct_reorder(LipidIon, MeanES), y = Exp_es))
a+
  # geom_point(aes(shape = Experiment, color = Cat), alpha = 0.8, size = 3) +
  geom_point(aes(shape = Experiment, color = Cat), alpha = 0.85, size = 3.5) +
  scale_shape_manual(values = c(13, 25, 18)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM), alpha = 0.75, width = 0.2, colour = "grey15",) +
  stat_summary(aes(x=LipidIon,y=MeanES), fun=mean, geom = "point", size=6, shape=20, alpha = 0.75, colour = "grey15",) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("seagreen", "tan2"))+
  coord_flip() + 
  labs(title = "DESI with in vivo and 2 out 3 invitro LC", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = all.ES.p.org %>% 
           filter(Sig == "Significant"), 
           aes(x = fct_reorder(LipidIon, MeanES), y = star.pos),
         pch=8, 
         size=1.5, stroke = 0.7, alpha = 0.75,
         colour="black") 
ggsave(filename = "./Figure_Panels/Fig.2d.pdf", width = 5, height = 5, useDingbats = FALSE)

