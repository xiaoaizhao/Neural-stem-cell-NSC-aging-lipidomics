### -------------------------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(ggthemes)
library("scales")
library(ggpubr)
library(ggbeeswarm)
library(stringi)
library(openxlsx)
library(eulerr)



## ----------------------------------------------------------------------------------------------------------------------------------
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Invitro_conc.lipid.qNSC.Rdata") #404 lipids
load("./Output_Data/GPMV_MsD.Raw.conc.Norm_Impt_291_lipids.Rdata") #291 lipids

Invitro.n <- MsD.lpd.rmv.abc(Invitro.q.conc)
GPMV <- GPMV.impt.norm.raw.conc.lpd.MsD %>% 
  rownames_to_column(var = "LipidIon")
GPMV.n <- MsD.lpd.rmv.abc(GPMV)



## ----------------------------------------------------------------------------------------------------------------------------------
Invitro.v.GPMV<- list(
  Invitro.lpd= unique(Invitro.n$LipidIon),
  GPMV.lpd= unique(GPMV.n$LipidIon))


## Venn diagram
pdf(file=paste0("./Figure_Panels/fig.S9C_left.pdf"), onefile=FALSE)
p1 <- euler(Invitro.v.GPMV)
plot(p1, fills = c("#d9e8e5",  "#E5E5E5"), 
     quantities = TRUE, legend = FALSE)
dev.off()


## ----------------------------------------------------------------------------------------------------------------------------------
common.GPMV.Invitro <- intersect(unique(Invitro.n$LipidIon), unique(GPMV.n$LipidIon))
save(common.GPMV.Invitro, file = "./Output_Data/Common.lipids_GPMV_Invitro.Rdata")


## ----------------------------------------------------------------------------------------------------------------------------------
rm(list = ls())
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/GPMV_MsD.Raw.conc.Norm_Impt_291_lipids.Rdata")
GPMV.conc <- GPMV.impt.norm.raw.conc.lpd.MsD %>% 
  rownames_to_column(var = "LipidIon") %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Concentration") 

load("./Output_Data/GPMV.MsD.TotalLipid.Rdata")

GPMV.lpd.pct <- left_join(GPMV.conc, GPMV.total.conc, by = "Sample") %>% 
  mutate(Lpd.pct = Concentration/TotalConc * 100)


## ----------------------------------------------------------------------------------------------------------------------------------
load("./Output_Data/Invitro_conc.lipid.qNSC.Rdata")
load("./Output_Data/Invitro.MsD.TotalLipid.A+Q.Rdata")

Invitro.conc <- MsD.lpd.rmv.abc(Invitro.q.conc) %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Concentration") %>% 
  filter(grepl("_qNSC-Q", Sample))

Invitro.q.total.conc <- InvitroAQ.total.conc %>% 
  filter(grepl("_qNSC-Q", Sample))

Invitro.lpd.pct <- left_join(Invitro.conc, Invitro.q.total.conc, by = "Sample") %>% 
  mutate(Lpd.pct = Concentration/TotalConc * 100)
  


## ----------------------------------------------------------------------------------------------------------------------------------
load("./Output_Data/Common.lipids_GPMV_Invitro.Rdata")

GPMV.cat <- MsD.lpd.rmv.abc(GPMV.lpd.pct) %>% 
  rowwise() %>% 
  mutate(Cat = ifelse(LipidIon %in% common.GPMV.Invitro, "Overlap", 
                      substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1))) %>% 
  mutate(Cat = ifelse(Cat =="HexCer", "Cer", Cat))

length(unique(GPMV.cat$LipidIon[GPMV.cat$Cat == "Overlap"])) #186

GPMV.sum.cat.sample <- GPMV.cat %>% 
  group_by(Cat, Sample) %>% 
  summarise(SumLpdPct = sum(Lpd.pct))

GPMV.mean.cat.pct.age <- GPMV.sum.cat.sample %>% 
  mutate(Age = ifelse(grepl("^O", Sample), "Old", "Young")) %>% 
  group_by(Age, Cat) %>% 
  summarise(MeanCatPct = mean(SumLpdPct)) %>% 
  mutate(Exp = "GPMV") %>% 
  mutate(Label = ifelse(Cat == "Overlap", Cat, paste0("Unique_", Cat)))
  


## ----------------------------------------------------------------------------------------------------------------------------------
Invitro.cat <- Invitro.lpd.pct %>% 
  rowwise() %>% 
  mutate(Cat = ifelse(LipidIon %in% common.GPMV.Invitro, "Overlap", 
                      substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)))

length(unique(Invitro.cat$LipidIon[Invitro.cat$Cat == "Overlap"])) #186

Invitro.sum.cat.sample <- Invitro.cat %>% 
  group_by(Cat, Sample) %>% 
  summarise(SumLpdPct = sum(Lpd.pct))

Invitro.mean.cat.pct.age <- Invitro.sum.cat.sample %>% 
  mutate(Age = ifelse(grepl("^O", Sample), "Old", "Young")) %>% 
  group_by(Age, Cat) %>% 
  summarise(MeanCatPct = mean(SumLpdPct)) %>% 
  mutate(Exp = "In vitro") %>% 
  mutate(Label = ifelse(Cat == "Overlap", Cat, paste0("Unique_", Cat)))


## ----------------------------------------------------------------------------------------------------------------------------------
Ovlp.lpd.pct <- bind_rows(GPMV.mean.cat.pct.age, Invitro.mean.cat.pct.age)

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")

mycolors <- lut.march25 %>% 
  filter(lipid_cat %in%Ovlp.lpd.pct$Cat) %>% 
  mutate(Label_cat = paste0("Unique_", lipid_cat))

Ovlp.lpd.pct$Label <- factor(Ovlp.lpd.pct$Label, levels = c("Overlap", mycolors$Label_cat))
Ovlp.lpd.pct <- Ovlp.lpd.pct %>% 
  arrange(Label)
Ovlp.lpd.pct$Age <- factor(Ovlp.lpd.pct$Age, levels = c("Young", "Old"))
Ovlp.lpd.pct$Exp <- factor(Ovlp.lpd.pct$Exp, levels = c("GPMV", "In vitro"))
b <- ggplot(Ovlp.lpd.pct, aes(x=Age, y=MeanCatPct))
b+geom_bar(aes(fill=Label),
           position="stack", stat="identity")+
  labs(title  = "GPMV vs. whole cell by molar concentration")+
  ylab("mol% of total lipid concentration")+
  scale_fill_manual(values = c("grey29", as.character(mycolors$Clr_list.1.18.))) +
  theme(legend.position= "right")+
  guides(fill=guide_legend(title="Cat")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black")) +
  facet_wrap(~Exp)
ggsave(filename = "./Figure_Panels/fig.S9C.middle.pdf", width = 5, height = 5, useDingbats = FALSE)


## ----------------------------------------------------------------------------------------------------------------------------------
rm(list = ls())
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Lpd.MsD.GPMV_Age_ES.Rdata")
load("./Output_Data/Lpd.MsD.Invitro_Qui_Age_ES.w.RT.Rdata")


## ----------------------------------------------------------------------------------------------------------------------------------
Invitro <- MsD.lpd.rmv.abc(Invitro.Q.lpd.ES.w.RT) %>% 
  ungroup() %>% 
  mutate(Cat = ifelse(es_g>0, "Up", "Down")) %>% 
  group_split(Cat)

Invitro.down <- Invitro[[1]] #184
Invitro.up <- Invitro[[2]] #220

## ----------------------------------------------------------------------------------------------------------------------------------
GPMV <- MsD.lpd.rmv.abc(GPMV.lpd.ES.w.RT) %>% 
  ungroup() %>% 
  mutate(Cat = ifelse(es_g>0, "Up", "Down")) %>% 
  group_split(Cat)

GPMV.down <- GPMV[[1]] #127
GPMV.up <- GPMV[[2]] #164


## ----------------------------------------------------------------------------------------------------------------------------------
GPMV.Invitro.up<- list(
  Invitro.lpd.up= unique(Invitro.up$LipidIon),
  GPMV.lpd.up= unique(GPMV.up$LipidIon))


## Venn diagram
pdf(file=paste0("./Figure_Panels/fig.S9C_right-middle.pdf"), onefile=FALSE)
p1 <- euler(GPMV.Invitro.up)
plot(p1, fills = c("#d9e8e5",  "#E5E5E5"), 
     quantities = TRUE, legend = FALSE) #78
dev.off()

overlap.up.invitro.GPMV <- intersect(unique(Invitro.up$LipidIon), unique(GPMV.up$LipidIon))
save(overlap.up.invitro.GPMV, file = "./Output_Data/Overlap_up.in.old.invitro.GPMV.Rdata")
## ----------------------------------------------------------------------------------------------------------------------------------
GPMV.Invitro.down<- list(
  Invitro.lpd.down= unique(Invitro.down$LipidIon),
  GPMV.lpd.down= unique(GPMV.down$LipidIon))


## Venn diagram
pdf(file=paste0("./Figure_Panels/fig.S9C_right.pdf"), onefile=FALSE)
p1 <- euler(GPMV.Invitro.down)
plot(p1, fills = c("#d9e8e5",  "#E5E5E5"), 
     quantities = TRUE, legend = FALSE) #46
dev.off()

overlap.down.invitro.GPMV <- intersect(unique(Invitro.down$LipidIon), unique(GPMV.down$LipidIon))
save(overlap.down.invitro.GPMV, file = "./Output_Data/Overlap_down.in.old.invitro.GPMV.Rdata")
