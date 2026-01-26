### Overlapping lipid analysis between In vitro and In vitro Experiment #2----------------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(ggthemes)
library("scales")
library(ggpubr)
library(ggbeeswarm)
library(stringi)
library(openxlsx)
library(eulerr)


## -------------------------------------------------------------------------------------------------------------------------
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Invitro_conc.lipid.qNSC.Rdata") #404 lipids
load("./Output_Data/Exp2_278_conc_lipid.MSD.format.Rdata") #278 lipids

Invitro.no.abc <- MsD.lpd.rmv.abc(Invitro.q.conc)

E2.no.abc <- MsD.lpd.rmv.abc(E2MsD.frmt)



## -------------------------------------------------------------------------------------------------------------------------
Invitro.v.Exp2<- list(
  Invitro.lpd= unique(Invitro.no.abc$LipidIon),
  Exp2.lpd= unique(E2.no.abc$LipidIon))


## Venn diagram
pdf(file=paste0("./Figure_Panels/EDFig.1b_left.pdf"), onefile=FALSE)
p1 <- euler(Invitro.v.Exp2)
plot(p1, fills = c("#d8f0f6", "#d9e8e5"), 
     quantities = TRUE, legend = FALSE)
dev.off()


## -------------------------------------------------------------------------------------------------------------------------
common.E2.E3 <- intersect(unique(Invitro.no.abc$LipidIon), unique(E2.no.abc$LipidIon))
save(common.E2.E3, file = "./Output_Data/Common.lipids_E2_E3.Rdata")


## -------------------------------------------------------------------------------------------------------------------------
rm(list = ls())
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Exp2_278_conc_lipid.MSD.format.Rdata")
E2.conc <- E2MsD.frmt %>% 
  select(LipidIon, contains("_N")) %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Concentration") 

load("./Output_Data/Exp2.KO.MsD.TotalLipid.Rdata")

E2.lpd.pct <- left_join(E2.conc, Exp2.total.conc, by = "Sample") %>% 
  mutate(Lpd.pct = Concentration/TotalConc * 100)


## -------------------------------------------------------------------------------------------------------------------------
load("./Output_Data/Invitro_conc.lipid.qNSC.Rdata")
load("./Output_Data/Invitro.QNSC.MsD.TotalLipid.No.Cholesterol.Rdata")

Invitro.conc <- Invitro.q.conc %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Concentration") %>% 
  filter(!grepl("ST", LipidIon))

Invitro.lpd.pct <- left_join(Invitro.conc, InvitroQ.total.conc.no.chol, by = "Sample") %>% 
  mutate(Lpd.pct = Concentration/TotalConc * 100)
  


## -------------------------------------------------------------------------------------------------------------------------
load("./Output_Data/Common.lipids_E2_E3.Rdata")

Exp2.cat <- MsD.lpd.rmv.abc(E2.lpd.pct) %>% 
  rowwise() %>% 
  mutate(Cat = ifelse(LipidIon %in% common.E2.E3, "Overlap", 
                      substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)))

length(unique(Exp2.cat$LipidIon[Exp2.cat$Cat == "Overlap"])) #175

Exp2.sum.cat.sample <- Exp2.cat %>% 
  group_by(Cat, Sample) %>% 
  summarise(SumLpdPct = sum(Lpd.pct))

Exp2.mean.cat.pct.age <- Exp2.sum.cat.sample %>% 
  mutate(Age = ifelse(grepl("^O", Sample), "Old", "Young")) %>% 
  group_by(Age, Cat) %>% 
  summarise(MeanCatPct = mean(SumLpdPct)) %>% 
  mutate(Exp = "In vitro Exp2") %>% 
  mutate(Label = ifelse(Cat == "Overlap", Cat, paste0("Unique_", Cat)))
  


## -------------------------------------------------------------------------------------------------------------------------
Invitro.cat <- MsD.lpd.rmv.abc(Invitro.lpd.pct) %>% 
  rowwise() %>% 
  mutate(Cat = ifelse(LipidIon %in% common.E2.E3, "Overlap", 
                      substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1))) %>% 
  filter(grepl("_qNSC-Q", Sample))

length(unique(Invitro.cat$LipidIon[Invitro.cat$Cat == "Overlap"])) #175

Invitro.sum.cat.sample <- Invitro.cat %>% 
  group_by(Cat, Sample) %>% 
  summarise(SumLpdPct = sum(Lpd.pct))

Invitro.mean.cat.pct.age <- Invitro.sum.cat.sample %>% 
  mutate(Age = ifelse(grepl("^O", Sample), "Old", "Young")) %>% 
  group_by(Age, Cat) %>% 
  summarise(MeanCatPct = mean(SumLpdPct)) %>% 
  mutate(Exp = "In vitro") %>% 
  mutate(Label = ifelse(Cat == "Overlap", Cat, paste0("Unique_", Cat)))


## -------------------------------------------------------------------------------------------------------------------------
Ovlp.lpd.pct <- bind_rows(Exp2.mean.cat.pct.age, Invitro.mean.cat.pct.age)

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")

mycolors <- lut.march25 %>% 
  filter(lipid_cat %in%Ovlp.lpd.pct$Cat) %>% 
  mutate(Label_cat = paste0("Unique_", lipid_cat))

Ovlp.lpd.pct$Label <- factor(Ovlp.lpd.pct$Label, levels = c("Overlap", mycolors$Label_cat))
Ovlp.lpd.pct <- Ovlp.lpd.pct %>% 
  arrange(Label)
Ovlp.lpd.pct$Age <- factor(Ovlp.lpd.pct$Age, levels = c("Young", "Old"))
Ovlp.lpd.pct$Exp <- factor(Ovlp.lpd.pct$Exp, levels = c("In vitro", "In vitro Exp2"))

b <- ggplot(Ovlp.lpd.pct, aes(x=Age, y=MeanCatPct))
b+geom_bar(aes(fill=Label),
           position="stack", stat="identity")+
  labs(title  = "2 invitro studies overlap by molar concentration")+
  ylab("mol% of total lipid concentration")+
  scale_fill_manual(values = c("grey29", as.character(mycolors$Clr_list.1.18.))) +
  theme(legend.position= "right")+
  guides(fill=guide_legend(title="Cat")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black")) +
  facet_wrap(~Exp)
ggsave(filename = "./Figure_Panels/EDFig.1b_middle.pdf", width = 5, height = 5, useDingbats = FALSE)


## -------------------------------------------------------------------------------------------------------------------------
rm(list = ls())
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Lpd.MsD.Age.ES_Exp2_all_KO.w.RT.Rdata")
load("./Output_Data/Lpd.MsD.Invitro_Qui_Age_ES.w.RT.Rdata")


## -------------------------------------------------------------------------------------------------------------------------
Invitro <- MsD.lpd.rmv.abc(Invitro.Q.lpd.ES.w.RT) %>% 
  ungroup() %>% 
  mutate(Cat = ifelse(es_g>0, "Up", "Down")) %>% 
  group_split(Cat)

Invitro.down <- Invitro[[1]] #184
Invitro.up <- Invitro[[2]] #220


## -------------------------------------------------------------------------------------------------------------------------
E2 <- MsD.lpd.rmv.abc(Exp2.lpd.es.g.MsD.allKO.wRT) %>% 
  filter(KO == "N") %>% 
  ungroup() %>% 
  mutate(Cat = ifelse(es_g>0, "Up", "Down")) %>% 
  group_split(Cat)

E2.down <- E2[[1]] #103
E2.up <- E2[[2]] #175


## -------------------------------------------------------------------------------------------------------------------------
Invitro.v.Exp2up<- list(
  Invitro.lpd.up= unique(Invitro.up$LipidIon),
  Exp2.lpd.up= unique(E2.up$LipidIon))


## Venn diagram
pdf(file=paste0("./Figure_Panels/EDFig.1b_right_middle.pdf"), onefile=FALSE)
p1 <- euler(Invitro.v.Exp2up)
plot(p1, fills = c("#d8f0f6", "#d9e8e5"), 
     quantities = TRUE, legend = FALSE) 
dev.off()

overlap.up.invitro <- intersect(unique(Invitro.up$LipidIon), unique(E2.up$LipidIon))
save(overlap.up.invitro, file = "./Output_Data/Overlap_up.in.old.2.invitro.studies.Rdata")
## -------------------------------------------------------------------------------------------------------------------------
Invitro.v.Exp2.down<- list(
  Invitro.lpd.down= unique(Invitro.down$LipidIon),
  Exp2.lpd.down= unique(E2.down$LipidIon))


## Venn diagram
pdf(file=paste0("./Figure_Panels/EDFig.1b_right.pdf"), onefile=FALSE)
p1 <- euler(Invitro.v.Exp2.down)
plot(p1, fills = c("#d8f0f6", "#d9e8e5"), 
     quantities = TRUE, legend = FALSE) 
dev.off()

overlap.dw.invitro <- intersect(unique(Invitro.down$LipidIon), unique(E2.down$LipidIon))
save(overlap.dw.invitro, file = "./Output_Data/Overlap_down.in.old.2.invitro.studies.Rdata")