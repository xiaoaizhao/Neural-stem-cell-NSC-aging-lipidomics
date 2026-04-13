
setwd(rstudioapi::getActiveProject())
library(tidyverse)
library(ggrepel)
rm(list = ls())
source("./Scripts/Function_scripts/Pre-processing_functions.R")
source("./Scripts/Function_scripts/Effect_size_functions.R")


load("./Output_Data/Exp2_SC.es.by.KO.in.old.compare.to.CTRL.Rdata")
load("./Output_Data/KO_LUT_paper.order2024-02-08.Rdata")
  
### ==== Side chains that show opposing effects in KO ====
### ==== old ====
# In KO dataset, positive effect size means higher in Mboat2 KO, negative effect size means higher in Control
MKO <- SC.all.ko.es.old %>% 
  filter(KO == "M")

### ==== lipids that show significant increase change in 2 of 2 invitro study + change in old is in the same direction as change in KO ====
load("./Output_Data/SC.2outof2LC.invitro.Padj.Sig.features.Rdata")
high.inold <- padj.sig.db. %>% 
  filter(MeanES>0)
high.in.old.ls <- unique(high.inold$Cla_SC) #51
high.inyoung <- padj.sig.db. %>% 
  filter(MeanES<0)
high.in.yng.ls <- unique(high.inyoung$Cla_SC) #10

KO.aging.sig <- inner_join(MKO, padj.sig.db., by = "Cla_SC") %>% 
  rename("es_g" = "es_g.x") %>% 
  filter(MeanES * es_g > 0) %>% 
  select(c(Cla_SC, es_g, se_g)) %>% 
  distinct() %>% 
  mutate(Aging.dir = case_when(
    Cla_SC %in% high.in.old.ls ~ "Higher in old",
    Cla_SC %in% high.in.yng.ls ~ "Higher in young",
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
  group_by(Cla_SC) %>% 
  filter(n_distinct(Sig) == 1) %>% 
  filter(Sig == "Significant")

padj.KO.SC.sig <- pval.from.CI.SC(CI.KO.sig, 0.05, 0) #7

padj.sig.KO.SC. <- padj.KO.SC.sig %>% #7
  filter(padj < 0.05) %>% 
  group_by(Cla_SC) %>% 
  group_modify(~{
    .x %>% 
      mutate(star.pos = case_when(
        es_g > 0 ~  max(es_g) + 2.4, 
        es_g < 0 ~ min(es_g) - 2)
      ) 
  }) %>% 
  arrange(es_g)
KO.responsive.SC.ls <- unique(padj.sig.KO.SC.$Cla_SC)

save(KO.responsive.SC.ls, file = "./Output_Data/KO_responsive_SC_list.Rdata")
### ==== Plot ====

up.p <-ggplot(padj.sig.KO.SC., aes(x = fct_reorder(Cla_SC, es_g), y = es_g))
up.p + 
  geom_point(shape = 11, colour = "#E15759", alpha = 0.85, size = 3) +
  geom_errorbar(aes(ymin = es_g - se_g, ymax = es_g + se_g, color = "black"), alpha = 1, width = 0.08) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 7.5), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "KO responsive SC (age significant)", x = "", y = "Effect size (Treatment vs. Control)") +
  geom_point(data = padj.sig.KO.SC. %>% 
               filter(Sig == "Significant"), 
             aes(x = Cla_SC, y = star.pos),
             pch=8, 
             size=1.5, stroke = 0.7, alpha = 0.75,
             colour="black") +
  coord_flip() 
ggsave("./Figure_Panels/EDFig.11e.PDF", width = 4, height = 4, useDingbats = FALSE)