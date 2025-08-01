
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)

source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_Data/Ef_Size_CONC.Lipid_InVivo.Rdata")
load("./Output_Data/Lipids.2outof3LC.invitro.Sig.features.Rdata")

Invivo <- ether.rename(Invivo.CONC.lpd.es.g)

CI95 <- function(df, mean, sem){
  df <- df %>% 
    rowwise() %>% 
    mutate(`Margin of Error` = 1.96 * {{sem}}) %>% 
    mutate(CI.lower = {{mean}} -`Margin of Error`) %>% 
    mutate(CI.upper = {{mean}} +`Margin of Error`)
}

c2.3.sig <- lpd.c2.in.3lc.df.sig %>% 
  select(LipidIon, MeanES) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
    .x %>% 
      summarise(SumES = unique(MeanES))
  })

Invivo.2of3.sig <- list(c2.3.sig, 
                  Invivo
                  ) %>% 
  reduce(inner_join, by = c("LipidIon")) %>%  #24 lipids overlap with In vivo
  rename("ES_summary" = "SumES") %>% 
  rename("ES_Invivo" = "es_g")


lpd.Invivo.2of3.sig.l <- Invivo.2of3.sig %>% 
  ungroup() %>% 
  select(LipidIon, starts_with("ES")) %>% 
  pivot_longer(-LipidIon, names_to = "Exp", values_to = "Effect size")

lpd.Invivo.2of3.sig <- lpd.Invivo.2of3.sig.l %>% 
  group_by(LipidIon) %>% 
  summarise(MeanES = mean(`Effect size`, na.rm = TRUE), SEM = sd(`Effect size`, na.rm = TRUE)/sqrt(sum(!is.na(`Effect size`)))) 
lpd.Invivo.2of3.sig.CI.LoUp <- CI95(lpd.Invivo.2of3.sig, MeanES, SEM)

c2of3sig.p <- inner_join(lpd.Invivo.2of3.sig.l, lpd.Invivo.2of3.sig.CI.LoUp, by = "LipidIon") %>% 
  mutate(Sig = ifelse( CI.lower > 0 & CI.upper > 0 | CI.lower < 0 & CI.upper < 0,
         "Significant", "Not significant")) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
    .x %>% 
      mutate(star.pos = case_when(
        Sig == "Significant" & MeanES > 0 ~  max(`Effect size`) + 0.8, 
        Sig == "Significant" & MeanES < 0 ~ min(`Effect size`) - 0.8)
)
  }) %>% 
  mutate(Exp = case_when(
    grepl("summary", Exp) ~ "In vitro summary",
    grepl("Invivo", Exp) ~ "In vivo",
  ))
  
c2of3sig.p$Exp <- factor(c2of3sig.p$Exp, levels = c( "In vitro summary", "In vivo"))

c2of3sig.p.sig <- c2of3sig.p %>% 
  filter(Sig == "Significant")

a <-ggplot(c2of3sig.p.sig, aes(x = fct_reorder(LipidIon, MeanES), y = `Effect size`))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 3.5) +
  scale_shape_manual(values = c(13, 25)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM), colour = "grey15", alpha = 0.75, width = 0.2) +
  stat_summary(aes(x=LipidIon,y=MeanES), fun=mean, geom = "point", size=3.5, shape=20, alpha = 0.75, colour = "grey15") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 8), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # scale_color_manual(values = c("cornflowerblue", "yellowgreen", "firebrick3"))+
  coord_flip() + 
  labs(title = "In vivo vs. 2 out of 3LC - significant lipids", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = c2of3sig.p.sig %>% 
           filter(Sig == "Significant"), 
           aes(x = fct_reorder(LipidIon, MeanES), y = star.pos),
         pch=8, 
         size=1.2, stroke = 0.7, alpha = 0.75,
         colour="black") +
  theme(legend.position = "bottom")
ggsave(filename = "./Figure_Panels/Fig.1j.pdf", width = 4, height = 5, useDingbats = FALSE)
