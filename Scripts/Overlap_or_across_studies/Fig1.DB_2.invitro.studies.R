
setwd(rstudioapi::getActiveProject())
rm(list = ls())

library(tidyverse)
library(ggpubr)
source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_Data/DB.MsD.Age.ES_Exp2_all_KO.Rdata") 
load("./Output_Data/DBPct.MsD.Invitro_Qui_Age_ES.Rdata") 

Exp2.ctrl.conc <- Exp2.MsD.DB.es.g.allKO %>% 
  filter(KO == "N") 

CI95 <- function(df, mean, sem){
  df <- df %>% 
    rowwise() %>% 
    mutate(`Margin of Error` = 1.96 * {{sem}}) %>% 
    mutate(CI.lower = {{mean}} -`Margin of Error`) %>% 
    mutate(CI.upper = {{mean}} +`Margin of Error`)
}

c2.in.2lc<- list(Invitro.Qui.MsD.CONC.DB.es.g, 
                   Exp2.ctrl.conc
                  ) %>% 
  reduce(full_join, by = c("Cla_DB")) #64

##replace NaN with NA
c2.NA <- c2.in.2lc %>% 
  select(Cla_DB, starts_with("es_g")) 

c2.NA.count <- c2.NA %>% 
  mutate(across(where(is.numeric), ~na_if(., NaN)))
  
c2.db <- c2.NA.count %>% 
  select(Cla_DB, starts_with("es_g")) %>% 
  mutate(Detected = 1 - rowSums(is.na(.))) %>% 
  filter(Detected >=1) #47 Cla_DB

c2.in.2lc.db <- c2.in.2lc %>% 
  filter(Cla_DB %in% c2.db$Cla_DB) %>% 
  select(Cla_DB, starts_with("es_g")) %>% 
  pivot_longer(-Cla_DB, names_to = "Exp", values_to = "Effect size") %>% 
  group_by(Cla_DB) %>% 
  summarise(MeanES = mean(`Effect size`, na.rm = TRUE), SEM = sd(`Effect size`, na.rm = TRUE)/sqrt(sum(!is.na(`Effect size`)))) 

c2.in.2lc.db.CI.LoUp <- CI95(c2.in.2lc.db, MeanES, SEM)

E2 <-Exp2.ctrl.conc %>% 
  filter(Cla_DB %in% c2.in.2lc.db.CI.LoUp$Cla_DB) %>% 
  select(Cla_DB, es_g) %>% 
  mutate(Exp = "In vitro Experiment #2")

Invitro <- Invitro.Qui.MsD.CONC.DB.es.g %>% 
  filter(Cla_DB %in% c2.in.2lc.db.CI.LoUp$Cla_DB) %>% 
  select(Cla_DB, es_g) %>% 
  mutate(Exp = "In vitro")


c2.in.2lc.df <- bind_rows(E2, Invitro)
p.c2.in.2lc.df <- left_join(c2.in.2lc.df, c2.in.2lc.db.CI.LoUp, by = "Cla_DB") %>% 
    mutate(Sat = case_when(
    as.numeric(substr(Cla_DB, nchar(Cla_DB), nchar(Cla_DB))) >1 ~ "PUFA",
    as.numeric(substr(Cla_DB, nchar(Cla_DB), nchar(Cla_DB))) ==1 ~ "MUFA",
    as.numeric(substr(Cla_DB, nchar(Cla_DB), nchar(Cla_DB))) ==0 ~ "SFA"
  )) %>% 
  mutate(Sig = ifelse( CI.lower > 0 & CI.upper > 0 | CI.lower < 0 & CI.upper < 0,
         "Significant", "Not significant"))
  

p.c2.in.2lc.df$Exp <- factor(p.c2.in.2lc.df$Exp, levels = c("In vitro", "In vitro Experiment #2"))
p.c2.in.2lc.df$Sat <- factor(p.c2.in.2lc.df$Sat, levels = c("SFA", "MUFA", "PUFA"))
db.c2.in.2lc <- p.c2.in.2lc.df

sig.db.c2.in.2lc <- db.c2.in.2lc %>% 
  filter(Sig == "Significant") #32 significant by CI95

nosig.db.c2.in.2lc <- db.c2.in.2lc %>% 
  filter(!Sig == "Significant") #62 not significant

## ==== Calculate p value based on CI ====
sig.db.c2.in.2lc.pval <- pval.from.CI.DB(sig.db.c2.in.2lc, 0.05, 0) #32

padj.sig.db. <- sig.db.c2.in.2lc.pval %>% 
  filter(padj < 0.05) %>% 
  group_by(Cla_DB) %>% 
  group_modify(~{
    .x %>% 
      mutate(star.pos = case_when(
        MeanES > 0 ~  max(es_g) + 0.8, 
        MeanES < 0 ~ min(es_g) - 0.8)
      ) 
  })

all.db.p <- bind_rows(nosig.db.c2.in.2lc, padj.sig.db.)
save(padj.sig.db., file = "./Output_Data/DB.2outof2LC.invitro.Padj.Sig.features.Rdata")


a <-ggplot(all.db.p, aes(x = fct_reorder(Cla_DB, MeanES), y = es_g))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 1.65, fill = "gold") +
  scale_shape_manual(values = c(17,2, 10)) +
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM, color = Sat), alpha = 0.85, width = 0.2) +
  stat_summary(aes(x=Cla_DB,y=MeanES, color = Sat), fun=mean, geom = "point", size=2, shape=20, alpha = 0.85) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 6), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("cornflowerblue", "yellowgreen", "firebrick3"))+
  coord_flip() +
  labs(title = "All DB 2 out of 2 LC in vitro", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = all.db.p %>%
           filter(Sig == "Significant"),
           aes(x = fct_reorder(Cla_DB, MeanES), y = star.pos),
         pch=8,
         size=0.9, stroke = 0.7, alpha = 0.75,
         colour="black") +
  theme(legend.position = "bottom")
ggsave(filename = "./Figure_Panels/EDFig.4b.pdf", width = 4, height = 6, useDingbats = FALSE)



a <-ggplot(padj.sig.db., aes(x = fct_reorder(Cla_DB, MeanES), y = es_g))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 2.5, fill = "gold") +
  scale_shape_manual(values = c(17,2, 10)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM, color = Sat), alpha = 0.95, width = 0.2, size = 1) +
  stat_summary(aes(x=Cla_DB,y=MeanES, color = Sat), fun=mean, geom = "point", size=3.5, shape=20, alpha = 0.85) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 6), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("cornflowerblue", "yellowgreen", "firebrick3"))+
  coord_flip() + 
  labs(title = "Padj Sig 2 out of 2 LC in vitro", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = padj.sig.db. %>% 
           filter(Sig == "Significant"), 
           aes(x = fct_reorder(Cla_DB, MeanES), y = star.pos),
         pch=8, 
         size=2, stroke = 0.7, alpha = 0.75,
         colour="black") 
ggsave(filename = "./Figure_Panels/Fig.1f.pdf", width = 5, height = 5, useDingbats = FALSE)
