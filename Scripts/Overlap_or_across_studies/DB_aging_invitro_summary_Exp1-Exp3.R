
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)

load("./Output_Data/Ef_Size_DB_pct_InVitro.Rdata") # Primary culture #1
load("./Output_Data/Ef_Size_CONC.DB_PCT_Exp2_all_KO.Rdata") # Primary culture #2
load("./Output_Data/Ef_Size_CONC.DB_pct_Exp3_Qui_aging.Rdata") # Primary culture #3

Exp2.ctrl.conc <- Exp2.CONC.DB.es.g.allKO %>% 
  filter(KO == "N") 

CI95 <- function(df, mean, sem){
  df <- df %>% 
    rowwise() %>% 
    mutate(`Margin of Error` = 1.96 * {{sem}}) %>% 
    mutate(CI.lower = {{mean}} -`Margin of Error`) %>% 
    mutate(CI.upper = {{mean}} +`Margin of Error`)
}

c2.in.3lc<- list(Exp3.Qui.CONC.DB.es.g, 
                   Exp2.ctrl.conc,
                  LC.Invitro.DB.es.g) %>% 
  reduce(full_join, by = c("Cla_DB")) #88
  
c2.db <- c2.in.3lc %>% 
  select(Cla_DB, starts_with("es_g")) %>% 
  mutate(Detected = 3 - rowSums(is.na(.))) %>% 
  filter(Detected >=2) #60 Cla_DB


c2.in.3lc.db <- c2.in.3lc %>% 
  filter(Cla_DB %in% c2.db$Cla_DB) %>% 
  select(Cla_DB, starts_with("es_g")) %>% 
  pivot_longer(-Cla_DB, names_to = "Exp", values_to = "Effect size") %>% 
  group_by(Cla_DB) %>% 
  summarise(MeanES = mean(`Effect size`, na.rm = TRUE), SEM = sd(`Effect size`, na.rm = TRUE)/sqrt(sum(!is.na(`Effect size`)))) 

c2.in.3lc.db.CI.LoUp <- CI95(c2.in.3lc.db, MeanES, SEM)

E1 <-LC.Invitro.DB.es.g %>% 
  filter(Cla_DB %in% c2.in.3lc.db.CI.LoUp$Cla_DB) %>% 
  select(Cla_DB, es_g) %>% 
  mutate(Exp = "Primary culture #1")

E2 <-Exp2.ctrl.conc %>% 
  filter(Cla_DB %in% c2.in.3lc.db.CI.LoUp$Cla_DB) %>% 
  select(Cla_DB, es_g) %>% 
  mutate(Exp = "Primary culture #2")

E3 <- Exp3.Qui.CONC.DB.es.g %>% 
  filter(Cla_DB %in% c2.in.3lc.db.CI.LoUp$Cla_DB) %>% 
  select(Cla_DB, es_g) %>% 
  mutate(Exp = "Primary culture #3")


c2.in.3lc.df <- bind_rows(E1, E2, E3)
p.c2.in.3lc.df <- left_join(c2.in.3lc.df, c2.in.3lc.db.CI.LoUp, by = "Cla_DB") %>% 
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
        Sig == "Significant" & MeanES > 0 ~  max(es_g) + 0.8, 
        Sig == "Significant" & MeanES < 0 ~ min(es_g) - 0.8)
)
  })
  

p.c2.in.3lc.df$Exp <- factor(p.c2.in.3lc.df$Exp, levels = c("Primary culture #1", "Primary culture #2", "Primary culture #3"))
p.c2.in.3lc.df$Sat <- factor(p.c2.in.3lc.df$Sat, levels = c("SFA", "MUFA", "PUFA"))
db.c2.in.3lc <- p.c2.in.3lc.df

save(db.c2.in.3lc, file = "./Output_Data/DB.2outof3LC.invitro.allfeatures.Rdata")

a <-ggplot(p.c2.in.3lc.df, aes(x = fct_reorder(Cla_DB, MeanES), y = es_g))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 1.65, fill = "gold") +
  scale_shape_manual(values = c(24, 2, 17, 10)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM, color = Sat), alpha = 0.85, width = 0.2) +
  stat_summary(aes(x=Cla_DB,y=MeanES, color = Sat), fun=mean, geom = "point", size=2, shape=20, alpha = 0.85) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 6), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("cornflowerblue", "yellowgreen", "firebrick3"))+
  coord_flip() + 
  labs(title = "2 out of all 3 LC in vitro", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = p.c2.in.3lc.df %>% 
           filter(Sig == "Significant"), 
           aes(x = fct_reorder(Cla_DB, MeanES), y = star.pos),
         pch=8, 
         size=0.9, stroke = 0.7, alpha = 0.75,
         colour="black") +
  theme(legend.position = "bottom")
ggsave(filename = "./Figure_Panels/EDFig.4b.pdf", width = 4, height = 6, useDingbats = FALSE)

p.c2.in.3lc.df.sig <- p.c2.in.3lc.df %>% 
  filter(Sig == "Significant")

save(p.c2.in.3lc.df.sig, file = "./Output_Data/DB.2outof3LC.invitro.Sig.features.Rdata")

a <-ggplot(p.c2.in.3lc.df.sig, aes(x = fct_reorder(Cla_DB, MeanES), y = es_g))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 2.5, fill = "gold") +
  scale_shape_manual(values = c(24, 2, 17, 10)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM, color = Sat), alpha = 0.95, width = 0.2, size = 1) +
  stat_summary(aes(x=Cla_DB,y=MeanES, color = Sat), fun=mean, geom = "point", size=3.5, shape=20, alpha = 0.85) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 6), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("cornflowerblue", "yellowgreen", "firebrick3"))+
  coord_flip() + 
  labs(title = "Significant 2 out of all 3 LC in vitro", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = p.c2.in.3lc.df.sig %>% 
           filter(Sig == "Significant"), 
           aes(x = fct_reorder(Cla_DB, MeanES), y = star.pos),
         pch=8, 
         size=2, stroke = 0.7, alpha = 0.75,
         colour="black") 
ggsave(filename = "./Figure_Panels/Fig.1g.pdf", width = 5, height = 5, useDingbats = FALSE)
