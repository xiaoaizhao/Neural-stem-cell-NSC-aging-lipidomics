### Class sum In vitro Experiment #2

setwd(rstudioapi::getActiveProject())
library(tidyverse)
rm(list=ls())

load("./Output_Data/Exp2_MsD.Norm_Impt_log2_conc_278_lipids.Rdata")

Exp2.class.sum <- 2^Exp2.log2.impt.norm.conc.MsD %>% 
  rownames_to_column(var = "LipidIon") %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Concentration") %>% 
  rowwise() %>% 
  mutate(Class = str_split(LipidIon, " ")[[1]][1]) %>% 
  group_by(Class, Sample) %>% 
  summarise(ClassSum = sum(Concentration))

save(Exp2.class.sum, file = "./Output_Data/Exp2.KO.MsD.ClassSum.Rdata")

Exp2.total.conc <- Exp2.class.sum %>% 
  group_by(Sample) %>% 
  summarise(TotalConc = sum(ClassSum))

save(Exp2.total.conc, file = "./Output_Data/Exp2.KO.MsD.TotalLipid.Rdata")
##====Bar plot====
rm(list = ls())
load("./Output_Data/Exp2.KO.MsD.ClassSum.Rdata")
load("./Output_Data/Exp2.KO.MsD.TotalLipid.Rdata")

### ==== control samples only ====
Exp2.cla.pct <- left_join(Exp2.class.sum, Exp2.total.conc, by = "Sample") %>% 
  filter(grepl("_N", Sample)) %>% 
  mutate(pct.per.cell = ClassSum / TotalConc *100) %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

Exp2.cla.pct$Age <- factor(Exp2.cla.pct$Age, levels = c("Young", "Old"))

Exp2.age.cla.pct <- Exp2.cla.pct %>% 
  group_by(Age, Class) %>% 
  summarise(Meanpct.per.cell = mean(pct.per.cell))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")

mycolors <- lut.march25 %>% 
  filter(lipid_cat %in%Exp2.age.cla.pct$Class)

Exp2.age.cla.pct$Class <- factor(Exp2.age.cla.pct$Class, levels = c(mycolors$lipid_cat))

Exp2.age.cla.pct <- Exp2.age.cla.pct %>% 
  arrange(Class)

b <- ggplot(Exp2.age.cla.pct, aes(x=Age, y=Meanpct.per.cell))
b+geom_bar(aes(fill=Class),
           position="stack", stat="identity")+
  labs(title  = "In vitro Exp2 - Lipid class composition")+
  ylab("mol% of total lipid concentration")+
  scale_fill_manual(values = as.character(mycolors$Clr_list.1.18.)) +
  theme(legend.position= "right")+
  guides(fill=guide_legend(title="Class")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"))
ggsave(filename = "./Figure_Panels/EDFig.1c.right.pdf", width = 5, height = 5, useDingbats = FALSE)
