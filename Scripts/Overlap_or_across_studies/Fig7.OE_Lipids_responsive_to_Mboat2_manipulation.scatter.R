## Results: only 2 significant lipids that are responsive to OE and also show significant age-related change.
setwd(rstudioapi::getActiveProject())
library(tidyverse)
library(ggrepel)
rm(list = ls())
source("./Scripts/Function_scripts/Pre-processing_functions.R")
source("./Scripts/Function_scripts/Effect_size_functions.R")

load("./Output_Data/M2OE_lpd.es.by.OE.in.Old.Rdata")

### ==== old ====
# In OE dataset, positive effect size means higher in Mboat2 OE, negative effect size means higher in Control
MOE <- M2.old.OE.ES

### ==== lipids that show significant increase change in 2 of 2 invitro study + change in old is in the same direction as change in KO ====
load("./Output_Data/Lipids.2of2LC.invitro.FDRSig.features.final.Rdata")
high.inold <- lpd.c2.in.2lc.df.sig.t.final %>% 
  filter(MeanES>0)
high.in.old.ls <- unique(high.inold$LipidIon)
high.inyoung <- lpd.c2.in.2lc.df.sig.t.final %>% 
  filter(MeanES<0)
high.in.yng.ls <- unique(high.inyoung$LipidIon)

OE.aging.sig <- inner_join(MOE, lpd.c2.in.2lc.df.sig.t.final, by = "LipidIon") %>% 
  rename("es_g" = "es_g.x") %>% 
  filter(MeanES * es_g < 0) %>% 
  select(c(LipidIon, es_g, se_g, MeanES, SEM)) %>% 
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
  group_by(LipidIon) %>% 
  filter(Sig == "Significant")


padj.OE.lpd.sig <- pval.from.CI(CI.OE.sig, 0.05, 0) %>% #2
  mutate(Aging.dir = case_when(
    LipidIon %in% high.in.old.ls ~ "Higher in old",
    LipidIon %in% high.in.yng.ls ~ "Higher in young",
  ))


OE.responsive.lpd.ls <- unique(padj.OE.lpd.sig$LipidIon)

# save(OE.responsive.lpd.ls, file = "./Output_Data/OE_responsive_lipid_list.Rdata")
# save(padj.OE.lpd.sig, file = "./Output_Data/OE_responsive.lpd.df.Rdata")
# 
# 
# load("./Output_Data/KO_responsive_lipid_list.Rdata")
# res.lpd.ls <- c(KO.responsive.lpd.ls, OE.responsive.lpd.ls)
# save(res.lpd.ls, file = "./Output_Data/Mboat2.responsive_lipid_list.Rdata")


padj.OE.lpd.sig$Aging.dir <- factor(padj.OE.lpd.sig$Aging.dir, levels = c("Higher in young", "Higher in old"))
### ==== Plot ====

pal4 <- c("cyan3", "magenta3")
a <- ggplot(padj.OE.lpd.sig, aes(es_g, MeanES))
a+
  geom_point(aes(color = Aging.dir), 
             size = 3, alpha = 0.85)+
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM), colour = "grey34", alpha = 0.4, width = 0.08) +
  geom_errorbar(aes(xmin = es_g - se_g, xmax = es_g + se_g), colour = "grey34", alpha = 0.4, width = 0.08) +
  scale_color_manual(values = pal4) +
  theme_classic()+
  labs(title = "OE responsive lipid (age significant)" , 
       y = "Effect size - (old vs. young)", 
       x = "Effect size - (OE vs. control)")+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(legend.position = "none") + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(aes(label = LipidIon), fontface = 'plain',
                  size = 3.5,colour = "black",
                  box.padding = unit(0.55, "lines"),
                  seed = 1234,
                  min.segment.length = 0,
                  max.overlaps = 45)
ggsave("./Figure_Panels/Fig.7C.PDF", width = 4, height = 4, useDingbats = FALSE)