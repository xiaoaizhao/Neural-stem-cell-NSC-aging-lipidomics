
setwd(rstudioapi::getActiveProject())
library(tidyverse)
library(ggrepel)
rm(list = ls())
source("./Scripts/Function_scripts/Pre-processing_functions.R")
source("./Scripts/Function_scripts/Effect_size_functions.R")

load("./Output_Data/Exp2_lpd.es.by.KO.in.old.compare.to.CTRL.Rdata")

### ==== old ====
# In KO dataset, positive effect size means higher in Mboat2 KO, negative effect size means higher in Control
MKO <- MsD.lpd.format(all.ko.es.old) %>% 
  filter(KO == "M")

### ==== lipids that show significant increase change in 2 of 2 invitro study + change in old is in the same direction as change in KO ====
load("./Output_Data/Lipids.2of2LC.invitro.FDRSig.features.final.Rdata")
high.inold <- lpd.c2.in.2lc.df.sig.t.final %>% 
  filter(MeanES>0)
high.in.old.ls <- unique(high.inold$LipidIon)
high.inyoung <- lpd.c2.in.2lc.df.sig.t.final %>% 
  filter(MeanES<0)
high.in.yng.ls <- unique(high.inyoung$LipidIon)

KO.aging.sig <- inner_join(MKO, lpd.c2.in.2lc.df.sig.t.final, by = "LipidIon") %>% 
  rename("es_g" = "es_g.x") %>% 
  filter(MeanES * es_g > 0) %>% 
  select(c(LipidIon, es_g, se_g)) %>% 
  distinct() %>% 
  mutate(Aging.dir = case_when(
    LipidIon %in% high.in.old.ls ~ "Higher in old",
    LipidIon %in% high.in.yng.ls ~ "Higher in young",
  ))


### ==== filter out to plot lipids with significant changes in either KO or OE condition ====
CI95 <- function(df, mean, sem){
  df <- df %>% 
    rowwise() %>% 
    mutate(`Margin of Error` = 1.96 * {{sem}}) %>% 
    mutate(CI.lower = {{mean}} -`Margin of Error`) %>% 
    mutate(CI.upper = {{mean}} +`Margin of Error`)
}

CI.KO.sig <- CI95(KO.aging.sig, es_g, se_g) %>% 
  mutate(Sig = ifelse( CI.lower > 0 & CI.upper > 0 | CI.lower < 0 & CI.upper < 0,
                       "Significant", "Not significant")) %>% 
  group_by(LipidIon) %>% 
  filter(Sig == "Significant")


padj.KO.lpd.sig <- pval.from.CI(CI.KO.sig, 0.05, 0) #9

padj.sig.KO.lpd. <- padj.KO.lpd.sig %>% #9
  filter(padj < 0.05) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
    .x %>% 
      mutate(star.pos = case_when(
        es_g > 0 ~  max(es_g) + 1.8, 
        es_g < 0 ~ min(es_g) - 1.8)
      ) 
  }) %>% 
  arrange(es_g)
KO.responsive.lpd.ls <- unique(padj.sig.KO.lpd.$LipidIon)

save(KO.responsive.lpd.ls, file = "./Output_Data/KO_responsive_lipid_list.Rdata")

save(padj.sig.KO.lpd., file = "./Output_Data/KO_responsive.lpd.df.Rdata")
### ==== Plot ====

up.p <-ggplot(padj.sig.KO.lpd., aes(x = fct_reorder(LipidIon, es_g), y = es_g))
up.p + 
  geom_point(shape = 11, colour = "#E15759", alpha = 0.85, size = 3) +
  geom_errorbar(aes(ymin = es_g - se_g, ymax = es_g + se_g, color = "black"), alpha = 1, width = 0.08) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 7.5), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "KO responsive lipid (age significant)", x = "", y = "Effect size (Treatment vs. Control)") +
  geom_point(data = padj.sig.KO.lpd. %>% 
               filter(Sig == "Significant"), 
             aes(x = LipidIon, y = star.pos),
             pch=8, 
             size=1.5, stroke = 0.7, alpha = 0.75,
             colour="black") +
  coord_flip() 
ggsave("./Figure_Panels/Fig.4f.PDF", width = 4, height = 4, useDingbats = FALSE)

