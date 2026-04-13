## Side chain analysis across 2 in vitro studies
## Remove LPC-O and LPE-O with 1 or 2 double bonds, as well as LPC-P and LPE-P with 0 or 1 double bond due to annotation ambiguity.

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)

source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_Data/SC.abundance.Invivo_Age_ES.Rdata")
load("./Output_Data/SC.abundance.Invitro_Qui_Age_ES.Rdata") 

CI95 <- function(df, mean, sem){
  df <- df %>% 
    rowwise() %>% 
    mutate(`Margin of Error` = 1.96 * {{sem}}) %>% 
    mutate(CI.lower = {{mean}} -`Margin of Error`) %>% 
    mutate(CI.upper = {{mean}} +`Margin of Error`)
}

Invivo <- Invivo.SC.es.g %>% 
  filter(!is.nan(es_g)) %>% 
  select(Cla_SC, es_g) %>% 
  rename("es_g.invivo" = "es_g")

Invitro <- Invitro.Qui.SC.es.g %>% 
  filter(!is.nan(es_g)) %>% 
  select(Cla_SC, es_g) %>% 
  rename("es_g.Exp3" = "es_g") %>% 
  filter(!(grepl("^LPC\\(O|^LPE\\(O", Cla_SC) & grepl(":1|:2", Cla_SC))) %>% 
  filter(!(grepl("^LPC\\(P|^LPE\\(P", Cla_SC) & grepl(":0|:1", Cla_SC))) #212

## only subset to include 
ovlp.db <- inner_join(Invitro, Invivo, by = "Cla_SC") %>% #46 overlapping DB features between in vivo and Exp3
  pivot_longer(-Cla_SC, names_to = "Exp", values_to = "Effect size")

Sum.es.db <- ovlp.db %>% 
  mutate(across(where(is.numeric), ~na_if(., NaN))) %>% 
  group_by(Cla_SC) %>% 
  summarise(MeanES = mean(`Effect size`, na.rm = TRUE), SEM = sd(`Effect size`, na.rm = TRUE)/sqrt(sum(!is.na(`Effect size`)))) 

DB.Invitro.Invivo.CI.LoUp <- CI95(Sum.es.db, MeanES, SEM)


c2of2all.p <- inner_join(ovlp.db, DB.Invitro.Invivo.CI.LoUp, by = "Cla_SC") %>% 
  mutate(Sat = ifelse(grepl("\\;", Cla_SC),
                      case_when(
                        as.numeric(substr(Cla_SC, str_locate(Cla_SC, "\\;")-1, str_locate(Cla_SC, "\\;")-1)) >1 ~ "PUFA",
                        as.numeric(substr(Cla_SC, str_locate(Cla_SC, "\\;")-1, str_locate(Cla_SC, "\\;")-1)) ==1 ~ "MUFA",
                        as.numeric(substr(Cla_SC, str_locate(Cla_SC, "\\;")-1, str_locate(Cla_SC, "\\;")-1)) ==0 ~ "SFA"
                      ),
                      case_when(
                        as.numeric(substr(Cla_SC, nchar(Cla_SC)-1, nchar(Cla_SC)-1)) >1 ~ "PUFA",
                        as.numeric(substr(Cla_SC, nchar(Cla_SC)-1, nchar(Cla_SC)-1)) ==1 ~ "MUFA",
                        as.numeric(substr(Cla_SC, nchar(Cla_SC)-1, nchar(Cla_SC)-1)) ==0 ~ "SFA"
                      ))) %>% 
  mutate(Sig = ifelse( CI.lower > 0 & CI.upper > 0 | CI.lower < 0 & CI.upper < 0,
                       "Significant", "Not significant")) %>% 
  group_by(Cla_SC) %>% 
  group_modify(~{
    .x %>% 
      mutate(star.pos = case_when(
        Sig == "Significant" & MeanES > 0 ~  max(`Effect size`) + 0.8, 
        Sig == "Significant" & MeanES < 0 ~ min(`Effect size`) - 0.8)
      )
  }) %>% 
  mutate(Exp = case_when(
    grepl("Exp3", Exp) ~ "In vitro",
    grepl("invivo", Exp) ~ "In vivo",
  ))


c2of2all.p$Exp <- factor(c2of2all.p$Exp, levels = c( "In vitro", "In vivo"))
c2of2all.p$Sat <- factor(c2of2all.p$Sat, levels = c("SFA", "MUFA", "PUFA"))


## ==== Calculate p value based on CI ====

a <-ggplot(c2of2all.p, aes(x = fct_reorder(Cla_SC, MeanES), y = `Effect size`))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.75, size = 4) +
  scale_shape_manual(values = c(17, 25)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM, color = Sat), alpha = 0.95, width = 0.3, size = 1.5) +
  stat_summary(aes(x=Cla_SC,y=MeanES), fun=mean, geom = "point", size=5, shape=20, alpha = 0.95, colour = "grey15") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 8), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("cornflowerblue", "yellowgreen", "firebrick3"))+
  coord_flip() + 
  labs(title = "In vitro vs. In vivo", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = c2of2all.p %>% 
               filter(Sig == "Significant"), 
             aes(x = fct_reorder(Cla_SC, MeanES), y = star.pos),
             pch=8, 
             size=2.5, stroke = 0.7, alpha = 0.75,
             colour="black") 
ggsave(filename = "./Figure_Panels/Fig.1i.pdf", width = 5, height = 5, useDingbats = FALSE)

## ==== highlight Mboat2 responsive Side chains====
load("./Output_Data/Mboat2.responsive_SideChain_list.Rdata")

invivo.ovlp.ls <- unique(c2of2all.p$Cla_SC)[unique(c2of2all.p$Cla_SC) %in% res.SC.ls]
invivo.ovlp.ls
# [1] "LPC(O-16:0)" "PE(20:4)"    "PE(O-18:1)" 
