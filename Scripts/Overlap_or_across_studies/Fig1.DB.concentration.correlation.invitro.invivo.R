
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)

load("./Output_Data/Invitro.Qui_CONC.DB_PCT_by_Class.Rdata") #Primary culture #3

Invitro.db <- MsD.Qui_CONC.DB.invitro %>% 
  mutate(Cla_DB = paste0(Class, DB_num)) %>% 
  group_by(Cla_DB, Age) %>% 
  summarise(mean.DBpct = mean(DB_Pct))

load("./Output_Data/MsD.Invivo_CONC.DB_PCT_by_Class.Rdata")

Invivo.db <- MsD.DB.pct.invivo %>% 
  mutate(Cla_DB = paste0(Class, DB_num)) %>% 
  group_by(Cla_DB, Age) %>% 
  summarise(mean.DBpct = mean(DB_Pct))

raw.DB.invitro.invivo <- inner_join(Invitro.db, Invivo.db, by = c("Age", "Cla_DB")) %>% 
  mutate(Class = substr(Cla_DB, 1, str_locate(Cla_DB, ":")-1)) %>% 
  rename("DBPct.Invitro" = "mean.DBpct.x") %>% 
  rename("DBPct.Invivo" = "mean.DBpct.y") 

raw.DB.invitro.invivo$Age <- factor(raw.DB.invitro.invivo$Age, levels = c("Young", "Old"))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")
mycolors <- lut.march25 %>% 
  filter(lipid_cat %in% raw.DB.invitro.invivo$Class)

mycolors$lipid_cat <- factor(mycolors$lipid_cat)

mycolors <- mycolors %>% 
  arrange(lipid_cat)

raw.DB.invitro.invivo$Class <- factor(raw.DB.invitro.invivo$Class, levels = levels(mycolors$lipid_cat))

a <- ggscatter(raw.DB.invitro.invivo, x = ("DBPct.Invitro"), y = ("DBPct.Invivo"), 
          color = "Class", shape = 16, size = 2.5,alpha = 0.9, # Points color, shape and size
          palette = as.character(mycolors$Clr_list.1.18.),
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.method = "pearson",
          xlab = "Double bond %mol Invitro", ylab = "Double bond %mol in vivo"
          )
a+ facet_wrap(~Age)+
  theme(text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))
ggsave(filename = "./Figure_Panels/EDFig.5e.pdf", width = 6, height = 6, useDingbats = FALSE)
