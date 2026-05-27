
setwd(rstudioapi::getActiveProject())
library(tidyverse)
library(ggrepel)
rm(list = ls())
source("./Scripts/Function_scripts/Pre-processing_functions.R")
source("./Scripts/Function_scripts/Effect_size_functions.R")

load("./Output_Data/M2OE_SC.es.by.OE.in.Old.Rdata")
  
### ==== Side chains that show opposing effects in OE ====
### ==== old ====
# In OE dataset, positive effect size means higher in Mboat2 oE, negative effect size means higher in Control
MOE <- M2.OE.old.SC.ES

### ==== lipids that show significant increase change in 2 of 2 invitro study + change in old is in the same direction as change in KO ====
load("./Output_Data/SC.2outof2LC.invitro.Padj.Sig.features.Rdata")
high.inold <- padj.sig.db. %>% 
  filter(MeanES>0)
high.in.old.ls <- unique(high.inold$Cla_SC) #51
high.inyoung <- padj.sig.db. %>% 
  filter(MeanES<0)
high.in.yng.ls <- unique(high.inyoung$Cla_SC) #10

OE.aging.sig <- inner_join(MOE, padj.sig.db., by = "Cla_SC") %>% 
  rename("es_g" = "es_g.x") %>% 
  filter(MeanES * es_g < 0) %>% 
  select(c(Cla_SC, es_g, se_g, MeanES, SEM)) %>% 
  distinct()


### ==== filter out to plot lipids with significant changes in either KO or OE condition ====
CI95 <- function(df, mean, sem){
  df <- df %>% 
    rowwise() %>% 
    mutate(`Margin of Error` = 1.96 * {{sem}}) %>% 
    mutate(CI.lower = {{mean}} -`Margin of Error`) %>% 
    mutate(CI.upper = {{mean}} +`Margin of Error`)
}

CI.OE.sig <- CI95(OE.aging.sig, es_g, se_g) %>% 
  mutate(Sig = ifelse( CI.lower > 0 & CI.upper > 0 | CI.lower < 0 & CI.upper < 0,
                       "Significant", "Not significant")) %>% 
  group_by(Cla_SC) %>% 
  filter(n_distinct(Sig) == 1) %>% 
  filter(Sig == "Significant")

padj.OE.SC.sig <- pval.from.CI.SC(CI.OE.sig, 0.05, 0) %>% #3
  mutate(Aging.dir = case_when(
    Cla_SC %in% high.in.old.ls ~ "Higher in old",
    Cla_SC %in% high.in.yng.ls ~ "Higher in young",
  ))

OE.responsive.SC.ls <- unique(padj.OE.SC.sig$Cla_SC)

save(OE.responsive.SC.ls, file = "./Output_Data/OE_responsive_SC_list.Rdata")

load("./Output_Data/KO_responsive_SC_list.Rdata")
res.SC.ls <- c(KO.responsive.SC.ls, OE.responsive.SC.ls)
save(res.SC.ls, file = "./Output_Data/Mboat2.responsive_SideChain_list.Rdata")


### ==== Plot ====
pal4 <- c("magenta3")
a <- ggplot(padj.OE.SC.sig, aes(es_g, MeanES))
a+
  geom_point(aes(color = Aging.dir), 
             size = 3, alpha = 0.85)+
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM), colour = "grey34", alpha = 0.4, width = 0.08) +
  geom_errorbar(aes(xmin = es_g - se_g, xmax = es_g + se_g), colour = "grey34", alpha = 0.4, width = 0.08) +
  scale_color_manual(values = pal4) +
  theme_classic()+
  labs(title = "OE responsive SC (age significant)" , 
       y = "Effect size - (old vs. young)", 
       x = "Effect size - (OE vs. control)")+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(aes(label = Cla_SC), fontface = 'plain',
                  size = 3.5,colour = "black",
                  box.padding = unit(0.55, "lines"),
                  seed = 1234,
                  min.segment.length = 0,
                  max.overlaps = 45)
ggsave("./Figure_Panels/fig.S12J.PDF", width = 4, height = 4, useDingbats = FALSE)