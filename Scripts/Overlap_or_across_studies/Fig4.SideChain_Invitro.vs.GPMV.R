## Side chain analysis across 2 in vitro studies
## Remove LPC-O and LPE-O with 1 or 2 double bonds, as well as LPC-P and LPE-P with 0 or 1 double bond due to annotation ambiguity.
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(ggrepel)
source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_Data/SC.abundance.Invitro_Qui_Age_ES.Rdata")
load("./Output_Data/SC.abundance.GPMV_Age_ES.Rdata")

CI95 <- function(df, mean, sem){
  df <- df %>% 
    rowwise() %>% 
    mutate(`Margin of Error` = 1.96 * {{sem}}) %>% 
    mutate(CI.lower = {{mean}} -`Margin of Error`) %>% 
    mutate(CI.upper = {{mean}} +`Margin of Error`)
}

GPMV <- GPMV.SC.es.g %>% 
  filter(!is.nan(es_g)) %>% 
  select(Cla_SC, es_g) %>% 
  rename("es_g.GPMV" = "es_g") %>% 
  filter(!grepl("^TG|^DG", Cla_SC)) %>%  #remove TG and DG classes from GPMV dataset
  filter(!(grepl("^LPC\\(O|^LPE\\(O", Cla_SC) & grepl(":1|:2", Cla_SC))) %>% 
  filter(!(grepl("^LPC\\(P|^LPE\\(P", Cla_SC) & grepl(":0|:1", Cla_SC)))

E3 <- Invitro.Qui.SC.es.g %>% 
  filter(!is.nan(es_g)) %>% 
  select(Cla_SC, es_g) %>% 
  rename("es_g.Exp3" = "es_g")

## only subset to include 
ovlp.db.E3GPMV <- inner_join(E3, GPMV, by = "Cla_SC") %>% 
  pivot_longer(-Cla_SC, names_to = "Exp", values_to = "Effect size")

Sum.es.db <- ovlp.db.E3GPMV %>% 
  mutate(across(where(is.numeric), ~na_if(., NaN))) %>% 
  group_by(Cla_SC) %>% 
  summarise(MeanES = mean(`Effect size`, na.rm = TRUE), SEM = sd(`Effect size`, na.rm = TRUE)/sqrt(sum(!is.na(`Effect size`)))) 

DB.GPMV.vs.Exp3.CI.LoUp <- CI95(Sum.es.db, MeanES, SEM)


c2of2all.p <- inner_join(ovlp.db.E3GPMV, DB.GPMV.vs.Exp3.CI.LoUp, by = "Cla_SC") %>% 
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
  mutate(Exp = case_when(
    grepl("Exp3", Exp) ~ "In vitro",
    grepl("GPMV", Exp) ~ "GPMV",
  ))


c2of2all.p$Exp <- factor(c2of2all.p$Exp, levels = c( "GPMV", "In vitro"))
c2of2all.p$Sat <- factor(c2of2all.p$Sat, levels = c("SFA", "MUFA", "PUFA"))


## ==== Calculate p value based on CI ====

sig.DB <- c2of2all.p %>% 
  filter(Sig == "Significant") #78
GPMV.Exp3.DB <- pval.from.CI.SC.summary(sig.DB, 0.05, 0) #78

GPMV.padj.sig <- GPMV.Exp3.DB %>% 
  filter(padj < 0.05) %>% #78
  group_by(Cla_SC) %>% 
  mutate(Cla_SC.l = ifelse(Sig == "Significant", 
                           paste0(Cla_SC, "<b>*</b>"), Cla_SC))

SC.sig.GPMV.invitro <- unique(GPMV.padj.sig$Cla_SC)
save(SC.sig.GPMV.invitro, file = "./Output_Data/FDR.sig.SC_GPMV_Invitro.Rdata")
### ==== get side chain feature of the top 30% ====
top30.df.hi <- GPMV.padj.sig %>%
  group_by(Cla_SC) %>%
  summarise(MeanES = unique(MeanES)) %>%
  mutate(percentile = percent_rank(MeanES)) %>%
  filter(percentile > 0.7) #12 side chain
top30.df.lo <- GPMV.padj.sig %>%
  filter(MeanES < 0) %>%
  group_by(Cla_SC) %>%
  summarise(MeanES = unique(MeanES)) %>%
  mutate(percentile = percent_rank(MeanES)) %>%
  filter(percentile < 0.3) #4 lipids
top30.ls <- c(top30.df.hi$Cla_SC, top30.df.lo$Cla_SC)
save(top30.ls, file = "./Output_Data/FDR.sig.SC_GPMV_Invitro.top30.pct.Rdata")
### ==== plot FDR significant ones for main figure ==== ####
a <-ggplot(GPMV.padj.sig, aes(x = fct_reorder(Cla_SC.l, MeanES), y = `Effect size`))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 3) +
  scale_shape_manual(values = c(15, 17)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM, color = Sat), alpha = 0.85, width = 0.3, size = 1) +
  stat_summary(aes(x=Cla_SC.l,y=MeanES), fun=mean, geom = "point", size=3, shape=20, alpha = 0.85, colour = "grey15") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 8), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("cornflowerblue", "yellowgreen", "firebrick3"))+
  coord_flip() + 
  labs(title = "SC GPMV vs. Exp3 FDR sig", x = "", y = "Effect size (Old vs. Young)") +
  theme(axis.text.y = element_markdown(colour = "black", size = 6))
ggsave(filename = paste0("./Figure_Panels/Fig.4C.pdf"), width = 5, height = 5, useDingbats = FALSE)


## ==== highlight Mboat2 responsive Side chains====
load("./Output_Data/Mboat2.responsive_SideChain_list.Rdata")

gpmv.ovlp.ls <- unique(GPMV.padj.sig$Cla_SC)[unique(GPMV.padj.sig$Cla_SC) %in% res.SC.ls]
gpmv.ovlp.ls
# [1] "PE(18:2)" "PE(20:4)" "PE(20:5)"


### ==== plot all overlapping DB for extended figure ==== ####
padj.ls <-GPMV.padj.sig %>% 
  select(Cla_SC, Exp, padj, Cla_SC.l)
SC.all <- c2of2all.p %>% 
  mutate(Cla_SC.l = ifelse(Sig == "Significant", 
                           paste0(Cla_SC, "<b>*</b>"), Cla_SC))

a <-ggplot(SC.all, aes(x = fct_reorder(Cla_SC.l, MeanES), y = `Effect size`))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 2) +
  scale_shape_manual(values = c(15, 17)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM, color = Sat), alpha = 0.85, width = 0.5, size = 0.5) +
  stat_summary(aes(x=Cla_SC.l,y=MeanES), fun=mean, geom = "point", size=2, shape=20, alpha = 0.85, colour = "grey15") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 6), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("cornflowerblue", "yellowgreen", "firebrick3"))+
  coord_flip() + 
  labs(title = "GPMV vs. Exp3 All SC", x = "", y = "Effect size (Old vs. Young)") +
  theme(axis.text.y = element_markdown(colour = "black", size = 6))+
  theme(legend.position = "bottom")
ggsave(filename = paste0("./Figure_Panels/fig.S9F.pdf"), width = 4, height = 16, useDingbats = FALSE)
