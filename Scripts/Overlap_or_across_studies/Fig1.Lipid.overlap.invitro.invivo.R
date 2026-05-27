
### ----------------------------------------------------------------------
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
source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Invitro_conc.lipid.qNSC.Rdata")
load("./Output_Data/Invivo.MsD.29.lipid.raw.conc.Rdata")

Invitro <- MsD.lpd.rmv.abc(Invitro.q.conc)
Invivo <- MsD.lpd.rmv.abc(Invivo.raw.conc.29.good.peak)


## -------------------------------------------------------------------------------------------------------------------------
Invitro.v.Invivo<- list(
  Invitro.lpd= unique(Invitro$LipidIon),
  Invivo.lpd= unique(Invivo$LipidIon))


## Venn diagram
pdf(file=paste0("./Figure_Panels/fig.S6B_left.pdf"), onefile=FALSE)
p1 <- euler(Invitro.v.Invivo)
plot(p1, fills = c("#d9e8e5", "#edd9ba"), 
     quantities = TRUE, legend = FALSE)
dev.off()


## -------------------------------------------------------------------------------------------------------------------------
common.Invitro.Invivo <- intersect(unique(Invivo$LipidIon), unique(Invitro$LipidIon))
save(common.Invitro.Invivo, file = "./Output_Data/Common.lipids_Invivo_Invitro.Rdata")
## -------------------------------------------------------------------------------------------------------------------------
rm(list = ls())

source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Invivo.MsD.29.lipid.raw.conc.Rdata")
Invivo.conc <- MsD.lpd.rmv.abc(Invivo.raw.conc.29.good.peak) %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Concentration") 

load("./Output_Data/Invivo.MsD.TotalLipid.Rdata")

Invivo.lpd.pct <- left_join(Invivo.conc, Invivo.total.conc, by = "Sample") %>% 
  mutate(Lpd.pct = Concentration/TotalConc * 100)


## -------------------------------------------------------------------------------------------------------------------------
load("./Output_Data/Invitro_conc.lipid.qNSC.Rdata")
load("./Output_Data/Invitro.MsD.TotalLipid.A+Q.Rdata")

Invitro.conc <- MsD.lpd.rmv.abc(Invitro.q.conc) %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Concentration") %>% 
  filter(grepl("_qNSC-Q", Sample))

Invitro.q.total.conc <- InvitroAQ.total.conc %>% 
  filter(grepl("_qNSC-Q", Sample))
  
Invitro.lpd.pct <- left_join(Invitro.conc, Invitro.q.total.conc, by = "Sample") %>% 
  mutate(Lpd.pct = Concentration/TotalConc * 100)
  


## -------------------------------------------------------------------------------------------------------------------------
load("./Output_Data/Common.lipids_Invivo_Invitro.Rdata")

Invivo.cat <- Invivo.lpd.pct %>% 
  rowwise() %>% 
  mutate(Cat = ifelse(LipidIon %in% common.Invitro.Invivo, "Overlap", 
                      substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)))

length(unique(Invivo.cat$LipidIon[Invivo.cat$Cat == "Overlap"])) #15

Invivo.sum.cat.sample <- Invivo.cat %>% 
  group_by(Cat, Sample) %>% 
  summarise(SumLpdPct = sum(Lpd.pct))

Invivo.mean.cat.pct.age <- Invivo.sum.cat.sample %>% 
  mutate(Age = ifelse(grepl("^O", Sample), "Old", "Young")) %>% 
  group_by(Age, Cat) %>% 
  summarise(MeanCatPct = mean(SumLpdPct)) %>% 
  mutate(Exp = "Invivo") %>% 
  mutate(Label = ifelse(Cat == "Overlap", Cat, paste0("Unique_", Cat)))
  


## -------------------------------------------------------------------------------------------------------------------------
Invitro.cat <- Invitro.lpd.pct %>% 
  rowwise() %>% 
  mutate(Cat = ifelse(LipidIon %in% common.Invitro.Invivo, "Overlap", 
                      substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)))

length(unique(Invitro.cat$LipidIon[Invitro.cat$Cat == "Overlap"])) #15

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
Ovlp.lpd.pct <- bind_rows(Invivo.mean.cat.pct.age, Invitro.mean.cat.pct.age)

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")

mycolors <- lut.march25 %>% 
  filter(lipid_cat %in%Ovlp.lpd.pct$Cat) %>% 
  mutate(Label_cat = paste0("Unique_", lipid_cat))

Ovlp.lpd.pct$Label <- factor(Ovlp.lpd.pct$Label, levels = c("Overlap", mycolors$Label_cat))
Ovlp.lpd.pct <- Ovlp.lpd.pct %>% 
  arrange(Label)
Ovlp.lpd.pct$Age <- factor(Ovlp.lpd.pct$Age, levels = c("Young", "Old"))

b <- ggplot(Ovlp.lpd.pct, aes(x=Age, y=MeanCatPct))
b+geom_bar(aes(fill=Label),
           position="stack", stat="identity")+
  labs(title  = "Invitro Invivo overlap by mol conc")+
  ylab("mol% of total lipid concentration")+
  scale_fill_manual(values = c("grey29", as.character(mycolors$Clr_list.1.18.))) +
  theme(legend.position= "right")+
  guides(fill=guide_legend(title="Cat")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black")) +
  facet_wrap(~Exp)
ggsave(filename = paste0("./Figure_Panels/fig.S6B.middle.pdf"), width = 5, height = 5, useDingbats = FALSE)


## -------------------------------------------------------------------------------------------------------------------------
rm(list = ls())
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Lpd.MsD.Invivo_Age_ES.w.RT.Rdata")
load("./Output_Data/Lpd.MsD.Invitro_Qui_Age_ES.w.RT.Rdata")


## -------------------------------------------------------------------------------------------------------------------------
Invitro <- MsD.lpd.rmv.abc(Invitro.Q.lpd.ES.w.RT) %>% 
  ungroup() %>% 
  mutate(Cat = ifelse(es_g>0, "Up", "Down")) %>% 
  group_split(Cat)

Invitro.down <- Invitro[[1]] #184
Invitro.up <- Invitro[[2]] #220


## -------------------------------------------------------------------------------------------------------------------------
Invivo <- MsD.lpd.rmv.abc(Invivo.lpd.ES.w.RT) %>% 
  ungroup() %>% 
  mutate(Cat = ifelse(es_g>0, "Up", "Down")) %>% 
  group_split(Cat)

Invivo.down <- Invivo[[1]] #13
Invivo.up <- Invivo[[2]] #16


## -------------------------------------------------------------------------------------------------------------------------
Invitro.v.Invivoup<- list(
  Invitro.lpd.up= unique(Invitro.up$LipidIon),
  Invivo.lpd.up= unique(Invivo.up$LipidIon))


## Venn diagram
pdf(file=paste0("./Figure_Panels/fig.S6B_right-middle.pdf"), onefile=FALSE)
p1 <- euler(Invitro.v.Invivoup)
plot(p1, fills = c("#d9e8e5", "#edd9ba"), 
     quantities = TRUE, legend = FALSE) #92
dev.off()
overlap.up.invitro.invivo <- intersect(unique(Invitro.up$LipidIon), unique(Invivo.up$LipidIon))
save(overlap.up.invitro.invivo, file = "./Output_Data/Overlap_up.in.old.invitro.invivo.Rdata")

## -------------------------------------------------------------------------------------------------------------------------
Invitro.v.InvivoDown<- list(
  Invitro.lpd.down= unique(Invitro.down$LipidIon),
  Invivo.lpd.down= unique(Invivo.down$LipidIon))


## Venn diagram
pdf(file=paste0("./Figure_Panels/fig.S6B_right.pdf"), onefile=FALSE)
p1 <- euler(Invitro.v.InvivoDown)
plot(p1, fills = c("#d9e8e5", "#edd9ba"), 
     quantities = TRUE, legend = FALSE) #37
dev.off()

overlap.down.invitro.invivo <- intersect(unique(Invitro.down$LipidIon), unique(Invivo.down$LipidIon))
save(overlap.down.invitro.invivo, file = "./Output_Data/Overlap_down.in.old.invitro.invivo.Rdata")
