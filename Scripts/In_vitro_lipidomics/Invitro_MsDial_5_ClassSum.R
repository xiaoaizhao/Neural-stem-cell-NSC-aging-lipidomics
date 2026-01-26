### Class sum In vitro 
setwd(rstudioapi::getActiveProject())
library(tidyverse)
rm(list=ls())

load("./Output_Data/Invitro.MsD.404.log2.norm.impt.no.O8aNSC.Rdata")

InvitroAQ.class.sum <- 2^Invitro.no.O8A.MsD %>% 
  rownames_to_column(var = "LipidIon") %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Concentration") %>% 
  rowwise() %>% 
  mutate(Class = str_split(LipidIon, " ")[[1]][1]) %>% 
  mutate(Class = ifelse(Class == "ST", "Cholesterol", Class)) %>% 
  group_by(Class, Sample) %>% 
  summarise(ClassSum = sum(Concentration))

save(InvitroAQ.class.sum, file = "./Output_Data/Invitro.MsD.ClassSum.A+Q.Rdata")

InvitroAQ.total.conc <- InvitroAQ.class.sum %>% 
  group_by(Sample) %>% 
  summarise(TotalConc = sum(ClassSum))

save(InvitroAQ.total.conc, file = "./Output_Data/Invitro.MsD.TotalLipid.A+Q.Rdata")
##====Bar plot====
rm(list = ls())
load("./Output_Data/Invitro.MsD.ClassSum.A+Q.Rdata")
load("./Output_Data/Invitro.MsD.TotalLipid.A+Q.Rdata")


Invitro.cla.pct <- left_join(InvitroAQ.class.sum, InvitroAQ.total.conc, by = "Sample") %>% 
  mutate(pct.per.cell = ClassSum / TotalConc *100) %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old")) %>% 
  mutate(CellType = ifelse(grepl("_qNSC-Q", Sample), "qNSC", "aNSC"))

Invitro.cla.pct$Age <- factor(Invitro.cla.pct$Age, levels = c("Young", "Old"))

Invitro.cell.age.cla.pct <- Invitro.cla.pct %>% 
  group_by(Age, Class, CellType) %>% 
  summarise(Meanpct.per.cell = mean(pct.per.cell))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")

mycolors <- lut.march25 %>% 
  filter(lipid_cat %in%Invitro.cell.age.cla.pct$Class)

Invitro.cell.age.cla.pct$Class <- factor(Invitro.cell.age.cla.pct$Class, levels = c(mycolors$lipid_cat))
Invitro.cell.age.cla.pct$CellType <- factor(Invitro.cell.age.cla.pct$CellType, levels = c("qNSC", "aNSC"))
Invitro.cell.age.cla.pct <- Invitro.cell.age.cla.pct %>% 
  arrange(Class)


b <- ggplot(Invitro.cell.age.cla.pct, aes(x=Age, y=Meanpct.per.cell))
b+geom_bar(aes(fill=Class),
           position="stack", stat="identity")+
  labs(title  = "In vitro all class composition")+
  ylab("mol% of total lipid concentration")+
  scale_fill_manual(values = as.character(mycolors$Clr_list.1.18.)) +
  theme(legend.position= "right")+
  guides(fill=guide_legend(title="Class")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black")) +
  facet_wrap(~CellType)
ggsave(filename = "./Figure_Panels/EDFig.1a.pdf", width = 5.5, height = 5, useDingbats=FALSE)


## ==== class pct after removing Cholesterol qNSC only==== ============================================================================
## ==== class pct after removing Cholesterol qNSC only==== ============================================================================
rm(list=ls())

load("./Output_Data/Invitro.MsD.404.log2.norm.impt.no.O8aNSC.Rdata")

InvitroQ.class.sum.no.chol <- 2^Invitro.no.O8A.MsD %>% 
  rownames_to_column(var = "LipidIon") %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Concentration") %>% 
  rowwise() %>% 
  mutate(Class = str_split(LipidIon, " ")[[1]][1]) %>% 
  filter(!Class == "ST") %>% 
  filter(grepl("qNSC-Q", Sample)) %>% 
  group_by(Class, Sample) %>% 
  summarise(ClassSum = sum(Concentration))

InvitroQ.total.conc.no.chol <- InvitroQ.class.sum.no.chol %>% 
  group_by(Sample) %>% 
  summarise(TotalConc = sum(ClassSum))

save(InvitroQ.total.conc.no.chol, file = "./Output_Data/Invitro.QNSC.MsD.TotalLipid.No.Cholesterol.Rdata")

Q.cla.pct.no.chol <- left_join(InvitroQ.class.sum.no.chol, InvitroQ.total.conc.no.chol, by = "Sample") %>% 
  mutate(pct.per.cell = ClassSum / TotalConc *100) %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

Q.cla.pct.no.chol$Age <- factor(Q.cla.pct.no.chol$Age, levels = c("Young", "Old"))

Q.cla.pct.no.chol <- Q.cla.pct.no.chol %>% 
  group_by(Age, Class) %>% 
  summarise(Meanpct.per.cell = mean(pct.per.cell))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")

mycolors <- lut.march25 %>% 
  filter(lipid_cat %in%Q.cla.pct.no.chol$Class)

Q.cla.pct.no.chol$Class <- factor(Q.cla.pct.no.chol$Class, levels = c(mycolors$lipid_cat))
Q.cla.pct.no.chol <- Q.cla.pct.no.chol %>% 
  arrange(Class)


b <- ggplot(Q.cla.pct.no.chol, aes(x=Age, y=Meanpct.per.cell))
b+geom_bar(aes(fill=Class),
           position="stack", stat="identity")+
  labs(title  = "In vitro qNSC no cholesterol class composition")+
  ylab("mol% of total lipid concentration")+
  scale_fill_manual(values = as.character(mycolors$Clr_list.1.18.)) +
  theme(legend.position= "right")+
  guides(fill=guide_legend(title="Class")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black")) 
ggsave(filename = "./Figure_Panels/EDFig.1c.left.pdf", width = 5, height = 5, useDingbats = FALSE)
