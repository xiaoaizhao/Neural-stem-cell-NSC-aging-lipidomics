## Plot Lipid Ontology Enrichment analysis result from Primary NSC culture #1 comparing old vs. young qNSC
## Plot Lipid Ontology Enrichment analysis result between Primary NSC culture #1 whole cell dataset vs. GPMV


### Primary NSC culture #1 
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

Exp1.LION <- read.csv("./Input_Data/Exp1_LION-cellular_component_entrichment.csv",
                     stringsAsFactors = F)

Exp1.LION <- Exp1.LION %>%
  mutate(., Cat = case_when(
    ES>0 & FDR.q.value<0.05 ~ "High in old",
    ES<0 & FDR.q.value<0.05 ~ "High in young",
    FDR.q.value>=0.05 ~ "Not Significant"
  )) %>%
  mutate(., Discription.Cap = str_to_sentence(Discription)) %>%
  mutate(., Discription.Cap = ifelse(grepl("Endoplasmic", Discription.Cap), "Endoplasmic reticulum", Discription.Cap))

Exp1.LION$Cat <- factor(Exp1.LION$Cat, levels = c("Not Significant", "High in old","High in young"))

palette1 <- c("grey60", "maroon", "darkgoldenrod")

a <- ggplot(Exp1.LION, aes(x=fct_reorder(Discription.Cap, ES), y=ES, fill=Cat))
a+ geom_bar(width = 0.88, stat = "identity") +
  coord_flip()+
  scale_fill_manual(values = palette1) +
  theme_classic()+
  labs(y="Enrichment Score (Old vs. Young)", x="", title = "Lipid Ontology Enrichment", fill = "")+
  theme(legend.position = "none")+
  theme(text=element_text(size = 15, face = "plain"),
        axis.text = element_text(colour = "black"))
ggsave(filename = "./Figure_Panels/Fig_3a.pdf", width = 7, height = 5, useDingbats=FALSE)


### Primary NSC culture #1 vs. GPMV

WCvPM<- read.csv("./Input_Data/WC_v_GPMV_LION-enrichment-job1.csv",
                     stringsAsFactors = F)

WCvPM <- WCvPM %>%
  mutate(., Cat = case_when(
    ES>0 & FDR.q.value<0.05 ~ "High in GPMV",
    ES<0 & FDR.q.value<0.05 ~ "High in Whole Cell",
    FDR.q.value>=0.05 ~ "Not Significant"
  )) %>%
  mutate(., Discription.Cap = str_to_sentence(Discription)) %>%
  mutate(., Discription.Cap = ifelse(grepl("Endoplasmic", Discription.Cap), "Endoplasmic reticulum", Discription.Cap))

WCvPM$Cat <- factor(WCvPM$Cat, levels = c("Not Significant", "High in Whole Cell","High in GPMV"))

palette1 <- c("grey60", "chartreuse3", "deepskyblue2")

a <- ggplot(WCvPM, aes(x=fct_reorder(Discription.Cap, ES), y=ES, fill=Cat))
a+ geom_bar(width = 0.88, stat = "identity") +
  coord_flip()+
  scale_fill_manual(values = palette1) +
  theme_classic()+
  labs(y="Enrichment Score (GPMV vs. Whole Cell)", x="", title = "Lipid Ontology Enrichment", fill = "")+
  theme(legend.position = "none")+
  theme(text=element_text(size = 15, face = "plain"),
        axis.text = element_text(colour = "black"))
ggsave(filename = "./Figure_Panels/Fig_S5b.pdf", width = 7, height = 5, useDingbats=FALSE)
