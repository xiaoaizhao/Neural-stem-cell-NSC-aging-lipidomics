
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)

load("./Output_Data/Ef_Size_CONC_DB_pct_Invivo.Rdata")
load('./Output_Data/DB.2outof3LC.invitro.Sig.features.Rdata')

CI95 <- function(df, mean, sem){
  df <- df %>% 
    rowwise() %>% 
    mutate(`Margin of Error` = 1.96 * {{sem}}) %>% 
    mutate(CI.lower = {{mean}} -`Margin of Error`) %>% 
    mutate(CI.upper = {{mean}} +`Margin of Error`)
}

c2.3.all <- p.c2.in.3lc.df.sig %>% 
  select(Cla_DB, MeanES) %>% 
  group_by(Cla_DB) %>% 
  group_modify(~{
    .x %>% 
      summarise(SumES = unique(MeanES))
  }) #19 Cla_DB

Invivo.2of3 <- list(c2.3.all, 
                  Invivo.CONC.DB.es.g
                  ) %>% 
  reduce(inner_join, by = c("Cla_DB")) %>%  #8 Cla_DB overlap with In vivo
  rename("ES_summary" = "SumES") %>% 
  rename("ES_Invivo" = "es_g") %>% 
  filter(!is.nan(ES_Invivo)) #7


DB.Invivo.2of3.l <- Invivo.2of3 %>% 
  ungroup() %>% 
  select(Cla_DB, starts_with("ES")) %>% 
  pivot_longer(-Cla_DB, names_to = "Exp", values_to = "Effect size")

DB.Invivo.2of3 <- DB.Invivo.2of3.l %>% 
  group_by(Cla_DB) %>% 
  summarise(MeanES = mean(`Effect size`, na.rm = TRUE), SEM = sd(`Effect size`, na.rm = TRUE)/sqrt(sum(!is.na(`Effect size`)))) 

DB.Invivo.2of3.CI.LoUp <- CI95(DB.Invivo.2of3, MeanES, SEM)


c2of3all.p <- inner_join(DB.Invivo.2of3.l, DB.Invivo.2of3.CI.LoUp, by = "Cla_DB") %>% 
  mutate(Sat = case_when(
    as.numeric(substr(Cla_DB, nchar(Cla_DB), nchar(Cla_DB))) >1 ~ "PUFA",
    as.numeric(substr(Cla_DB, nchar(Cla_DB), nchar(Cla_DB))) ==1 ~ "MUFA",
    as.numeric(substr(Cla_DB, nchar(Cla_DB), nchar(Cla_DB))) ==0 ~ "SFA"
  )) %>% 
  mutate(Sig = ifelse( CI.lower > 0 & CI.upper > 0 | CI.lower < 0 & CI.upper < 0,
         "Significant", "Not significant")) %>% 
  group_by(Cla_DB) %>% 
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
  

c2of3all.p$Exp <- factor(c2of3all.p$Exp, levels = c( "In vitro summary", "In vivo"))
c2of3all.p$Sat <- factor(c2of3all.p$Sat, levels = c("SFA", "MUFA", "PUFA"))


a <-ggplot(c2of3all.p, aes(x = fct_reorder(Cla_DB, MeanES), y = `Effect size`))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.75, size = 4) +
  scale_shape_manual(values = c(13, 15)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM, color = Sat), alpha = 0.95, width = 0.3, size = 1.5) +
  stat_summary(aes(x=Cla_DB,y=MeanES), fun=mean, geom = "point", size=5, shape=20, alpha = 0.95, colour = "grey15") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 8), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("cornflowerblue", "yellowgreen", "firebrick3"))+
  coord_flip() + 
  labs(title = "In vivo vs. 2 out of 3LC - significant DB", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = c2of3all.p %>% 
           filter(Sig == "Significant"), 
           aes(x = fct_reorder(Cla_DB, MeanES), y = star.pos),
         pch=8, 
         size=2.5, stroke = 0.7, alpha = 0.75,
         colour="black") 
ggsave(filename = "./Figure_Panels/Fig.1l.pdf", width = 5, height = 5, useDingbats = FALSE)