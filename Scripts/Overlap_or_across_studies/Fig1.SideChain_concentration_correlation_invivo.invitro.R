
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)

load("./Output_Data/Invitro.SC.analysis.Rdata") 

Invitro.SC <- Invitro.SumSC %>% 
  mutate(Cla_SC = paste0(Class, "(", Sep.SC, ")")) %>% 
  mutate(Age = ifelse(grepl("^O", Samples), "Old", "Young")) %>% 
  group_by(Cla_SC, Age) %>% 
  summarise(meanSumSC = mean(SumSC))

load("./Output_Data/Invivo.SC.analysis.Rdata")

Invivo.SC <- Invivo.SC.sum %>% 
  mutate(Cla_SC = paste0(Class, "(", Sep.SC, ")")) %>% 
  mutate(Age = ifelse(grepl("^O", Sample), "Old", "Young")) %>% 
  group_by(Cla_SC, Age) %>% 
  summarise(meanSumSC = mean(SumSC))

raw.SC.invivo.invitro <- inner_join(Invitro.SC, Invivo.SC, by = c("Age", "Cla_SC")) %>% 
  mutate(Class = substr(Cla_SC, 1, str_locate(Cla_SC, "\\(")-1)) %>% 
  rename("SCsum.Invitro" = "meanSumSC.x") %>% 
  rename("SCsum.Invivo" = "meanSumSC.y") %>% 
  mutate(Log2Invitro = log2(SCsum.Invitro)) %>% 
  mutate(Log2Invivo = log2(SCsum.Invivo))

raw.SC.invivo.invitro$Age <- factor(raw.SC.invivo.invitro$Age, levels = c("Young", "Old"))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")
mycolors <- lut.march25%>% 
  filter(lipid_cat %in% raw.SC.invivo.invitro$Class)

mycolors$lipid_cat <- factor(mycolors$lipid_cat)

mycolors <- mycolors %>% 
  arrange(lipid_cat)

raw.SC.invivo.invitro$Class <- factor(raw.SC.invivo.invitro$Class, levels = levels(mycolors$lipid_cat))

a <- ggscatter(raw.SC.invivo.invitro, x = "SCsum.Invitro", y = "SCsum.Invivo", 
          color = "Class", shape = 16, size = 2.5,alpha = 0.9, # Points color, shape and size
          palette = as.character(mycolors$Clr_list.1.18.),
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.method = "pearson",
          xlab = "Side chain concentration (uM) In vitro", ylab = "Side chain concentration (uM) In vivo"
          )
a+ facet_wrap(~Age)+
  theme(text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))
ggsave(filename = "./Figure_Panels/fig.S6E.pdf", width = 6, height = 6, useDingbats = FALSE)
