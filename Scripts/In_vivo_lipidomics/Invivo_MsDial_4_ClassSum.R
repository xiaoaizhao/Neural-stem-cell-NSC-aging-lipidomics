### Class sum In vivo
setwd(rstudioapi::getActiveProject())
library(tidyverse)
rm(list=ls())

load("./Output_Data/Invivo_MsD.Norm_Impt_log2_conc_29.lipids.Rdata")

Invivo.raw.conc.29.good.peak <- 2^Invivo.log2.impt.norm.conc.MsD.good.peak %>% 
  rownames_to_column(var = "LipidIon")

save(Invivo.raw.conc.29.good.peak, file = "./Output_Data/Invivo.MsD.29.lipid.raw.conc.Rdata")

Invivo.class.sum <- Invivo.raw.conc.29.good.peak %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Concentration") %>% 
  rowwise() %>% 
  mutate(Class = str_split(LipidIon, " ")[[1]][1]) %>% 
  mutate(Class = ifelse(Class == "ST", "Cholesterol", Class)) %>% 
  group_by(Class, Sample) %>% 
  summarise(ClassSum = sum(Concentration))

save(Invivo.class.sum, file = "./Output_Data/Invivo.MsD.ClassSum.Rdata")

Invivo.total.conc <- Invivo.class.sum %>% 
  group_by(Sample) %>% 
  summarise(TotalConc = sum(ClassSum))

save(Invivo.total.conc, file = "./Output_Data/Invivo.MsD.TotalLipid.Rdata")
##====Bar plot====
rm(list = ls())
load("./Output_Data/Invivo.MsD.ClassSum.Rdata")
load("./Output_Data/Invivo.MsD.TotalLipid.Rdata")


Invivo.cla.pct<- left_join(Invivo.class.sum, Invivo.total.conc, by = "Sample") %>% 
  mutate(pct.per.cell = ClassSum / TotalConc *100) %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

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
  labs(title  = "In vivo all class composition")+
  ylab("mol% of total lipid concentration")+
  scale_fill_manual(values = as.character(mycolors$Clr_list.1.18.)) +
  theme(legend.position= "right")+
  guides(fill=guide_legend(title="Class")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black")) 
ggsave(filename = "./Figure_Panels/fig.S6A.pdf", width = 5, height = 5, useDingbats = FALSE)


