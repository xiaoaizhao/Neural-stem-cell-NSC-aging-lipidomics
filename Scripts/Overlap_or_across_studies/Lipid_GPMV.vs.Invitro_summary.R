

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_Data/Ef_Size_Conc.Lipid_GPMV.Rdata")
load("./Output_Data/Lipids.2outof3LC.invitro.Sig.features.Rdata")

GPMV <- ether.rename(GPMV.conc.lpd.es.g)

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

GPMV.2of3.sig <- list(c2.3.sig, 
                  GPMV
                  ) %>% 
  reduce(inner_join, by = c("LipidIon")) %>%  #91 lipids overlap with In vivo
  rename("ES_summary" = "SumES") %>% 
  rename("ES_GPMV" = "es_g")


lpd.GPMV.2of3.sig.l <- GPMV.2of3.sig %>% 
  ungroup() %>% 
  select(LipidIon, starts_with("ES")) %>% 
  pivot_longer(-LipidIon, names_to = "Exp", values_to = "Effect size")

lpd.GPMV.2of3.sig <- lpd.GPMV.2of3.sig.l %>% 
  group_by(LipidIon) %>% 
  summarise(MeanES = mean(`Effect size`, na.rm = TRUE), SEM = sd(`Effect size`, na.rm = TRUE)/sqrt(sum(!is.na(`Effect size`)))) 

lpd.GPMV.2of3.sig.CI.LoUp <- CI95(lpd.GPMV.2of3.sig, MeanES, SEM)


c2of3sig.p <- inner_join(lpd.GPMV.2of3.sig.l, lpd.GPMV.2of3.sig.CI.LoUp, by = "LipidIon") %>% 
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
    grepl("GPMV", Exp) ~ "GPMV",
  ))
  

c2of3sig.p$Exp <- factor(c2of3sig.p$Exp, levels = c( "In vitro summary", "GPMV"))

##==== All overlapping lipids between in vitro summary and GPMV ====
a <-ggplot(c2of3sig.p, aes(x = fct_reorder(LipidIon, MeanES), y = `Effect size`))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 2.5) +
  scale_shape_manual(values = c(13, 15)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM), colour = "grey15", alpha = 0.75, width = 0.2) +
  stat_summary(aes(x=LipidIon,y=MeanES), fun=mean, geom = "point", size=3.5, shape=20, alpha = 0.75, colour = "grey15") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 8), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() + 
  labs(title = "GPMV vs. 2 out of 3LC - significant lipids", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = c2of3sig.p %>% 
           filter(Sig == "Significant"), 
           aes(x = fct_reorder(LipidIon, MeanES), y = star.pos),
         pch=8, 
         size=1.2, stroke = 0.7, alpha = 0.75,
         colour="black") +
  theme(legend.position = "bottom")
ggsave(filename = "./Figure_Panels/EDFig.8f.pdf", width = 4, height = 11, useDingbats = FALSE)


##==== Take only top 50% of significant lipids between in vitro summary and GPMV for main figure ====
c2of3sig.p.sig <- c2of3sig.p %>% 
  filter(Sig == "Significant")


GPMV.invitro.top50 <- c2of3sig.p.sig %>% 
  group_by(LipidIon) %>% 
  summarise(MeanES = unique(MeanES)) %>% 
  mutate(percentile = percent_rank(abs(MeanES))) %>% 
  filter(percentile > 0.5) # 20 lipids

GPMV.invitro.top50.p <- c2of3sig.p.sig %>% 
  filter(LipidIon %in% GPMV.invitro.top50$LipidIon)

a <-ggplot(GPMV.invitro.top50.p, aes(x = fct_reorder(LipidIon, MeanES), y = `Effect size`))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 2.5) +
  scale_shape_manual(values = c(13, 15)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM), colour = "grey15", alpha = 0.75, width = 0.2) +
  stat_summary(aes(x=LipidIon,y=MeanES), fun=mean, geom = "point", size=3.5, shape=20, alpha = 0.75, colour = "grey15") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 8), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # scale_color_manual(values = c("cornflowerblue", "yellowgreen", "firebrick3"))+
  coord_flip() + 
  labs(title = "Top 50% SIG GPMV vs. 2 out of 3LC - significant lipids", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = GPMV.invitro.top50.p %>% 
           filter(Sig == "Significant"), 
           aes(x = fct_reorder(LipidIon, MeanES), y = star.pos),
         pch=8, 
         size=1.2, stroke = 0.7, alpha = 0.75,
         colour="black") +
  theme(legend.position = "bottom")
ggsave(filename = "./Figure_Panels/Fig.3c.pdf", width = 4, height = 5, useDingbats = FALSE)
