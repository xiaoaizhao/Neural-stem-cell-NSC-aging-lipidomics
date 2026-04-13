## Side chain analysis across 2 in vitro studies
## Remove LPC-O and LPE-O with 1 or 2 double bonds, as well as LPC-P and LPE-P with 0 or 1 double bond due to annotation ambiguity.
## Need to run the following 2 scripts first
### 'Fig4.KO_SideChain_responsive_to_Mboat2_manipulation.R'
### 'Fig4.OE_SideChain_responsive_to_Mboat2_manipulation.R'
setwd(rstudioapi::getActiveProject())
rm(list = ls())

library(tidyverse)
library(ggpubr)
source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_Data/SC.abundance.MsD.Age.ES_Exp2_all_KO.Rdata") 
load("./Output_Data/SC.abundance.Invitro_Qui_Age_ES.Rdata") 

Exp2.ctrl.conc <- Exp2.SC.es.g.allKO %>% 
  filter(KO == "N") %>% #145
  filter(!(grepl("^LPC\\(O|^LPE\\(O", Cla_SC) & grepl(":1|:2", Cla_SC))) %>% 
  filter(!(grepl("^LPC\\(P|^LPE\\(P", Cla_SC) & grepl(":0|:1", Cla_SC))) #137

CI95 <- function(df, mean, sem){
  df <- df %>% 
    rowwise() %>% 
    mutate(`Margin of Error` = 1.96 * {{sem}}) %>% 
    mutate(CI.lower = {{mean}} -`Margin of Error`) %>% 
    mutate(CI.upper = {{mean}} +`Margin of Error`)
}

c2.in.2lc<- list(Invitro.Qui.SC.es.g, 
                   Exp2.ctrl.conc
                  ) %>% 
  reduce(full_join, by = c("Cla_SC")) #234

##replace NaN with NA
c2.NA <- c2.in.2lc %>% 
  select(Cla_SC, starts_with("es_g")) 

c2.NA.count <- c2.NA %>% 
  mutate(across(where(is.numeric), ~na_if(., NaN)))
  
c2.db <- c2.NA.count %>% 
  select(Cla_SC, starts_with("es_g")) %>% 
  mutate(Detected = 1 - rowSums(is.na(.))) %>% 
  filter(Detected >=1) #120 Cla_DB

c2.in.2lc.db <- c2.in.2lc %>% 
  filter(Cla_SC %in% c2.db$Cla_SC) %>% 
  select(Cla_SC, starts_with("es_g")) %>% 
  pivot_longer(-Cla_SC, names_to = "Exp", values_to = "Effect size") %>% 
  group_by(Cla_SC) %>% 
  summarise(MeanES = mean(`Effect size`, na.rm = TRUE), SEM = sd(`Effect size`, na.rm = TRUE)/sqrt(sum(!is.na(`Effect size`)))) 

c2.in.2lc.db.CI.LoUp <- CI95(c2.in.2lc.db, MeanES, SEM)

E2 <-Exp2.ctrl.conc %>% 
  filter(Cla_SC %in% c2.in.2lc.db.CI.LoUp$Cla_SC) %>% 
  select(Cla_SC, es_g) %>% 
  mutate(Exp = "In vitro Experiment #2")

Invitro <- Invitro.Qui.SC.es.g %>% 
  filter(Cla_SC %in% c2.in.2lc.db.CI.LoUp$Cla_SC) %>% 
  select(Cla_SC, es_g) %>% 
  mutate(Exp = "In vitro")


c2.in.2lc.df <- bind_rows(E2, Invitro)
p.c2.in.2lc.df <- left_join(c2.in.2lc.df, c2.in.2lc.db.CI.LoUp, by = "Cla_SC") %>% 
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
         "Significant", "Not significant"))
  

p.c2.in.2lc.df$Exp <- factor(p.c2.in.2lc.df$Exp, levels = c("In vitro", "In vitro Experiment #2"))
p.c2.in.2lc.df$Sat <- factor(p.c2.in.2lc.df$Sat, levels = c("SFA", "MUFA", "PUFA"))
db.c2.in.2lc <- p.c2.in.2lc.df

sig.db.c2.in.2lc <- db.c2.in.2lc %>% 
  filter(Sig == "Significant") #116 significant by CI95

nosig.db.c2.in.2lc <- db.c2.in.2lc %>% 
  filter(!Sig == "Significant") #124 not significant

## ==== Calculate p value based on CI ====
sig.db.c2.in.2lc.pval <- pval.from.CI.SC(sig.db.c2.in.2lc, 0.05, 0) #116

padj.sig.db. <- sig.db.c2.in.2lc.pval %>% 
  filter(padj < 0.05) %>% 
  group_by(Cla_SC) %>% 
  group_modify(~{
    .x %>% 
      mutate(star.pos = case_when(
        MeanES > 0 ~  max(es_g) + 0.8, 
        MeanES < 0 ~ min(es_g) - 0.8)
      ) 
  })

all.db.p <- bind_rows(nosig.db.c2.in.2lc, padj.sig.db.)
save(padj.sig.db., file = "./Output_Data/SC.2outof2LC.invitro.Padj.Sig.features.Rdata")


a <-ggplot(all.db.p, aes(x = fct_reorder(Cla_SC, MeanES), y = es_g))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 1.65, fill = "gold") +
  scale_shape_manual(values = c(17,2, 10)) +
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM, color = Sat), alpha = 0.85, width = 0.2) +
  stat_summary(aes(x=Cla_SC,y=MeanES, color = Sat), fun=mean, geom = "point", size=2, shape=20, alpha = 0.85) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 6), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("cornflowerblue", "yellowgreen", "firebrick3"))+
  coord_flip() +
  labs(title = "All DB 2 out of 2 LC in vitro", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = all.db.p %>%
           filter(Sig == "Significant"),
           aes(x = fct_reorder(Cla_SC, MeanES), y = star.pos),
         pch=8,
         size=0.9, stroke = 0.7, alpha = 0.75,
         colour="black") +
  theme(legend.position = "bottom")
ggsave(filename = "./Figure_Panels/EDFig.4b.pdf", width = 4, height = 10, useDingbats = FALSE)




#### ==== All side chain features in the top 30% ====
top30.all.hi <- all.db.p %>% 
  group_by(Cla_SC) %>% 
  summarise(MeanES = unique(MeanES)) %>% 
  mutate(percentile = percent_rank(MeanES)) %>% 
  filter(percentile > 0.7) #38 side chain
top30.all.lo <- all.db.p %>% 
  filter(MeanES < 0) %>% 
  group_by(Cla_SC) %>% 
  summarise(MeanES = unique(MeanES)) %>% 
  mutate(percentile = percent_rank(MeanES)) %>% 
  filter(percentile < 0.3) #11 lipids
top30.all.ls <- c(top30.all.hi$Cla_SC, top30.all.lo$Cla_SC)

top30.all.sc <- all.db.p %>% 
  filter(Cla_SC %in% top30.all.ls)

a <-ggplot(top30.all.sc, aes(x = fct_reorder(Cla_SC, MeanES), y = es_g))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 1.65, fill = "gold") +
  scale_shape_manual(values = c(17,2, 10)) +
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM, color = Sat), alpha = 0.85, width = 0.2) +
  stat_summary(aes(x=Cla_SC,y=MeanES, color = Sat), fun=mean, geom = "point", size=2, shape=20, alpha = 0.85) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 6), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("cornflowerblue", "yellowgreen", "firebrick3"))+
  coord_flip() +
  labs(title = "Top 30 percentile All SC 2 in vitro", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = top30.all.sc %>%
               filter(Sig == "Significant"),
             aes(x = fct_reorder(Cla_SC, MeanES), y = star.pos),
             pch=8,
             size=0.7, stroke = 0.7, alpha = 0.75,
             colour="black") +
  theme(legend.position = "bottom")
ggsave(filename = "./Figure_Panels/Fig.1f.pdf", width = 4, height = 8, useDingbats = FALSE)


## ==== highlight Mboat2 responsive Side chains====
load("./Output_Data/Mboat2.responsive_SideChain_list.Rdata")

invitro.ovlp.ls <- unique(top30.all.sc$Cla_SC)[unique(top30.all.sc$Cla_SC) %in% res.SC.ls]
invitro.ovlp.ls
# [1] "LPC(O-16:0)" "PE(20:4)"    "PE(20:5)"    "PE(O-18:1)"  "PI(20:5)"  