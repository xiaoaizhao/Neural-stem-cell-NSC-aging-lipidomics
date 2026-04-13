
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

load("./Output_Data/E2.SC.analysis.Rdata")

Exp2.ctrl.SC <- E2.SumSC %>% 
  mutate(KO = substr(Samples, nchar(Samples)-1, nchar(Samples))) %>% 
  filter(KO == "_N") %>% 
  mutate(Cla_SC = paste0(Class, "(", Sep.SC, ")")) %>% 
  mutate(Age = ifelse(grepl("^O", Samples), "Old", "Young")) %>% 
  group_by(Cla_SC, Age) %>% 
  summarise(meanSumSC = mean(SumSC))

raw.SC.2.invitro <- inner_join(Invitro.SC, Exp2.ctrl.SC, by = c("Age", "Cla_SC")) %>% 
  mutate(Class = substr(Cla_SC, 1, str_locate(Cla_SC, "\\(")-1)) %>% 
  rename("SCsum.Invitro" = "meanSumSC.x") %>% 
  rename("SCsum.InvitroExp2" = "meanSumSC.y") %>% 
  mutate(Log2Invitro = log2(SCsum.Invitro)) %>% 
  mutate(Log2Invitro_Exp2 = log2(SCsum.InvitroExp2))

raw.SC.2.invitro$Age <- factor(raw.SC.2.invitro$Age, levels = c("Young", "Old"))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")
mycolors <- lut.march25%>% 
  filter(lipid_cat %in% raw.SC.2.invitro$Class)

mycolors$lipid_cat <- factor(mycolors$lipid_cat)

mycolors <- mycolors %>% 
  arrange(lipid_cat)

raw.SC.2.invitro$Class <- factor(raw.SC.2.invitro$Class, levels = levels(mycolors$lipid_cat))

a <- ggscatter(raw.SC.2.invitro, x = "Log2Invitro", y = "Log2Invitro_Exp2", 
          color = "Class", shape = 16, size = 2.5,alpha = 0.9, # Points color, shape and size
          palette = as.character(mycolors$Clr_list.1.18.),
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.method = "pearson",
          xlab = "Log2 side chain concentration In vitro", ylab = "Log2 side chain concentration In vitro Experiment #2"
          )
a+ facet_wrap(~Age)+
  theme(text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))
ggsave(filename = "./Figure_Panels/EDFig.4a.pdf", width = 6, height = 6, useDingbats = FALSE)
