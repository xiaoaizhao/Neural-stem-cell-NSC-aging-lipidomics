
setwd(rstudioapi::getActiveProject())
library(tidyverse)
rm(list = ls())

##====Primary culture #3====
load("./Output_Data/M2PM.CTRL.ClassConcSum.Rdata")

M2PM.total.conc <- M2PM.ctrl.cla.sum %>% 
  group_by(Samples) %>% 
  summarise(TotalConc = sum(ClassSum))

M2PM.cla.pct <- left_join(M2PM.ctrl.cla.sum, M2PM.total.conc, by = "Samples") %>% 
  mutate(pct.per.cell = ClassSum / TotalConc *100)

M2PM.cla.pct$Age <- factor(M2PM.cla.pct$Age, levels = c("Young", "Old"))

M2PM.age.cla.pct <- M2PM.cla.pct %>% 
  group_by(Age, Class) %>% 
  summarise(Meanpct.per.cell = mean(pct.per.cell))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")

mycolors <- lut.march25 %>% 
  filter(lipid_cat %in%M2PM.age.cla.pct$Class)

M2PM.age.cla.pct$Class <- factor(M2PM.age.cla.pct$Class, levels = c(mycolors$lipid_cat))

M2PM.age.cla.pct <- M2PM.age.cla.pct %>% 
    arrange(Class)

b <- ggplot(M2PM.age.cla.pct, aes(x=Age, y=Meanpct.per.cell))
  b+geom_bar(aes(fill=Class),
  position="stack", stat="identity")+
    labs(title  = "Primary cultuer #3")+
    ylab("mol% of total lipid concentration")+
    scale_fill_manual(values = as.character(mycolors$Clr_list.1.18.)) +
    theme(legend.position= "right")+
    guides(fill=guide_legend(title="Class")) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"))
  ggsave(filename = "./Figure_Panels/EDFig.1a.pdf", width = 5, height = 5, useDingbats = FALSE)


##====Lipidyzer====
rm(list = ls())
load("./Output_Data/Lipidyzer_qNSC.ClassConcSum.dupTG.removed.Rdata")

Ldz.noFFA <- Cla.ldz.Q.nodup %>% 
  filter(!Class == "FFA")

Ldz.total.conc <- Ldz.noFFA %>% 
  group_by(Samples) %>% 
  summarise(TotalConc = sum(ClassSum))

Ldz.cla.pct <- left_join(Ldz.noFFA, Ldz.total.conc, by = "Samples") %>% 
  mutate(pct.per.cell = ClassSum / TotalConc *100) %>% 
  mutate(Class = ifelse(Class == "CE", "ChE", Class))

ldz.age.cla.pct <- Ldz.cla.pct %>% 
  group_by(Age, Class) %>% 
  summarise(Meanpct.per.cell = mean(pct.per.cell))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")

mycolors <- lut.march25 %>% 
  filter(lipid_cat %in%ldz.age.cla.pct$Class)
mycolors$lipid_cat <- factor(mycolors$lipid_cat)

ldz.age.cla.pct$Class <- factor(ldz.age.cla.pct$Class, levels = c(mycolors$lipid_cat))

ldz.age.cla.pct <- ldz.age.cla.pct %>% 
    arrange(Class)

ldz.age.cla.pct$Age <- factor(ldz.age.cla.pct$Age, levels = c("Young", "Old"))

b <- ggplot(ldz.age.cla.pct, aes(x=Age, y=Meanpct.per.cell))
  b+geom_bar(aes(fill=Class),
  position="stack", stat="identity")+
    labs(title  = "Lipid class composition - Lipidyzer")+
    ylab("mol% of total lipid concentration")+
    scale_fill_manual(values = as.character(mycolors$Clr_list.1.18.)) +
    theme(legend.position= "right")+
    guides(fill=guide_legend(title="Class")) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"))
    ggsave(filename = "./Figure_Panels/EDFig.3b.pdf", width = 5, height = 5, useDingbats = FALSE)

##====In vivo====
load("./Output_Data/Invivo.ClassConcSum.Rdata")   
Invivo.total.conc <- Invivo.Cla.org.l %>% 
  group_by(Samples) %>% 
  summarise(TotalConc = sum(ClassSum))

Invivo.cla.pct <- left_join(Invivo.Cla.org.l, Invivo.total.conc, by = "Samples") %>% 
  mutate(pct.per.cell = ClassSum / TotalConc *100)

Invivo.cla.pct$Age <- factor(Invivo.cla.pct$Age, levels = c("Young", "Old"))

Invivo.cla.pct <- Invivo.cla.pct %>% 
  group_by(Age, Class) %>% 
  summarise(Meanpct.per.cell = mean(pct.per.cell))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")

mycolors <- lut.march25 %>% 
  filter(lipid_cat %in%Invivo.cla.pct$Class)

Invivo.cla.pct$Class <- factor(Invivo.cla.pct$Class, levels = c(mycolors$lipid_cat))

Invivo.cla.pct <- Invivo.cla.pct %>% 
  arrange(Class)

b <- ggplot(Invivo.cla.pct, aes(x=Age, y=Meanpct.per.cell))
b+geom_bar(aes(fill=Class),
           position="stack", stat="identity")+
  labs(title  = "In vivo")+
  ylab("mol% of total lipid concentration")+
  scale_fill_manual(values = as.character(mycolors$Clr_list.1.18.)) +
  theme(legend.position= "right")+
  guides(fill=guide_legend(title="Class")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"))
ggsave(filename = "./Figure_Panels/EDFig.5a.pdf", width = 5, height = 5, useDingbats = FALSE)

##====GPMV====
rm(list = ls())
load("./Output_Data/GPMV_ClassSum_CONC.483_lipids.Rdata")

GPMV.total.conc <- GPMV.conc.classsum %>% 
  group_by(Sample) %>% 
  summarise(TotalConc = sum(Class_sum))

G.cla.pct <- left_join(GPMV.conc.classsum, GPMV.total.conc, by = "Sample") %>% 
  mutate(pct.per.cell = Class_sum / TotalConc *100) %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

G.cla.pct$Age <- factor(G.cla.pct$Age, levels = c("Young", "Old"))

GPMV.age.cla.pct <- G.cla.pct %>% 
  group_by(Age, Class) %>% 
  summarise(Meanpct.per.cell = mean(pct.per.cell))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")

mycolors <- lut.march25 %>% 
  filter(lipid_cat %in%GPMV.age.cla.pct$Class)

GPMV.age.cla.pct$Class <- factor(GPMV.age.cla.pct$Class, levels = c(mycolors$lipid_cat))

GPMV.age.cla.pct <- GPMV.age.cla.pct %>% 
    arrange(Class)

b <- ggplot(GPMV.age.cla.pct, aes(x=Age, y=Meanpct.per.cell))
  b+geom_bar(aes(fill=Class),
  position="stack", stat="identity")+
    labs(title  = "Lipid class composition - GPMV")+
    ylab("mol% of total lipid concentration")+
    scale_fill_manual(values = as.character(mycolors$Clr_list.1.18.)) +
    theme(legend.position= "right")+
    guides(fill=guide_legend(title="Class")) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"))
  ggsave(filename = "./Figure_Panels/EDFig.8b.pdf", width = 5, height = 5, useDingbats = FALSE)
