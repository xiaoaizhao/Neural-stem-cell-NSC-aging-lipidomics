
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)

load("./Output_Data/Invitro.Qui_CONC.DB_PCT_by_Class.Rdata") 

Invitro.db <- MsD.Qui_CONC.DB.invitro %>% 
  mutate(Cla_DB = paste0(Class, DB_num)) %>% 
  group_by(Cla_DB, Age) %>% 
  summarise(mean.DBpct = mean(DB_Pct))

load("./Output_Data/MsD.Exp2.KO_CONC.DB_PCT_by_Class.Rdata")

Exp2.ctrl.db <- MsD.E2_CONC.DB %>% 
  mutate(KO = substr(Sample, nchar(Sample)-1, nchar(Sample))) %>% 
  filter(KO == "_N") %>% 
  mutate(Cla_DB = paste0(Class, DB_num)) %>% 
  group_by(Cla_DB, Age) %>% 
  summarise(mean.DBpct = mean(DB_Pct))

raw.DB.2.invitro <- inner_join(Invitro.db, Exp2.ctrl.db, by = c("Age", "Cla_DB")) %>% 
  mutate(Class = substr(Cla_DB, 1, str_locate(Cla_DB, ":")-1)) %>% 
  rename("DBPct.Invitro" = "mean.DBpct.x") %>% 
  rename("DBPct.InvitroExp2" = "mean.DBpct.y") 

raw.DB.2.invitro$Age <- factor(raw.DB.2.invitro$Age, levels = c("Young", "Old"))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")
mycolors <- lut.march25%>% 
  filter(lipid_cat %in% raw.DB.2.invitro$Class)

mycolors$lipid_cat <- factor(mycolors$lipid_cat)

mycolors <- mycolors %>% 
  arrange(lipid_cat)

raw.DB.2.invitro$Class <- factor(raw.DB.2.invitro$Class, levels = levels(mycolors$lipid_cat))

a <- ggscatter(raw.DB.2.invitro, x = ("DBPct.Invitro"), y = ("DBPct.InvitroExp2"), 
          color = "Class", shape = 16, size = 2.5,alpha = 0.9, # Points color, shape and size
          palette = as.character(mycolors$Clr_list.1.18.),
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.method = "pearson",
          xlab = "Double bond %mol In vitro", ylab = "Double bond %mol In vitro Experiment #2"
          )
a+ facet_wrap(~Age)+
  theme(text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))
ggsave(filename = "./Figure_Panels/EDFig.4a.pdf", width = 6, height = 6, useDingbats = FALSE)
