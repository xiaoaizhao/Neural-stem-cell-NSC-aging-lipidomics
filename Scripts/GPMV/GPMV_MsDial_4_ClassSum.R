### Class sum GPMV
setwd(rstudioapi::getActiveProject())
library(tidyverse)
rm(list=ls())

load("./Output_Data/GPMV_MsD.Raw.conc.Norm_Impt_291_lipids.Rdata")

GPMV.class.sum <- GPMV.impt.norm.raw.conc.lpd.MsD %>% 
  rownames_to_column(var = "LipidIon") %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Concentration") %>% 
  rowwise() %>% 
  mutate(Class = str_split(LipidIon, " ")[[1]][1]) %>% 
  mutate(Class = ifelse(Class == "ST", "Cholesterol", Class)) %>% 
  mutate(Class = ifelse(Class == "HexCer", "Cer", Class)) %>% 
  group_by(Class, Sample) %>% 
  summarise(ClassSum = sum(Concentration))

save(GPMV.class.sum, file = "./Output_Data/GPMV.MsD.ClassSum.Rdata")

GPMV.total.conc <- GPMV.class.sum %>% 
  group_by(Sample) %>% 
  summarise(TotalConc = sum(ClassSum))

save(GPMV.total.conc, file = "./Output_Data/GPMV.MsD.TotalLipid.Rdata")
##====Bar plot====
rm(list = ls())
load("./Output_Data/GPMV.MsD.ClassSum.Rdata")
load("./Output_Data/GPMV.MsD.TotalLipid.Rdata")


GPMV.cla.pct <- left_join(GPMV.class.sum, GPMV.total.conc, by = "Sample") %>% 
  mutate(pct.per.cell = ClassSum / TotalConc *100) %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

GPMV.cla.pct$Age <- factor(GPMV.cla.pct$Age, levels = c("Young", "Old"))

GPMV.age.cla.pct <- GPMV.cla.pct %>% 
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
  labs(title  = "GPMV all class composition")+
  ylab("mol% of total lipid concentration")+
  scale_fill_manual(values = as.character(mycolors$Clr_list.1.18.)) +
  theme(legend.position= "right")+
  guides(fill=guide_legend(title="Class")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black")) 
ggsave(filename = "./Figure_Panels/fig.S9B.pdf", width = 5, height = 5, useDingbats = FALSE)


