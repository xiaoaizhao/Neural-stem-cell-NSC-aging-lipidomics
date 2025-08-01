
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)

load("./Output_Data/Exp3.Qui_CONC.DB_PCT_by_Class.Rdata") #Primary culture #3

E3.db <- Qui_CONC.DB.E3 %>% 
  mutate(Cla_DB = paste0(Class, DB_num)) %>% 
  group_by(Cla_DB, Age) %>% 
  summarise(mean.DBpct = mean(DB_Pct))

load("./Output_Data/Invivo_CONC_DB_PCT_by_Class.Rdata")

Invivo.db <- Invivo.CONC.DB.PCT %>% 
  mutate(Cla_DB = paste0(Class, DB_num)) %>% 
  group_by(Cla_DB, Age) %>% 
  summarise(mean.DBpct = mean(DB_Pct))

raw.DB.E3.invivo <- inner_join(E3.db, Invivo.db, by = c("Age", "Cla_DB")) %>% 
  mutate(Class = substr(Cla_DB, 1, str_locate(Cla_DB, ":")-1)) %>% 
  rename("DBPct.Exp3" = "mean.DBpct.x") %>% 
  rename("DBPct.Invivo" = "mean.DBpct.y") 

raw.DB.E3.invivo$Age <- factor(raw.DB.E3.invivo$Age, levels = c("Young", "Old"))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")
mycolors <- lut.march25 %>% 
  filter(lipid_cat %in% raw.DB.E3.invivo$Class)

mycolors$lipid_cat <- factor(mycolors$lipid_cat)

mycolors <- mycolors %>% 
  arrange(lipid_cat)

raw.DB.E3.invivo$Class <- factor(raw.DB.E3.invivo$Class, levels = levels(mycolors$lipid_cat))

a <- ggscatter(raw.DB.E3.invivo, x = ("DBPct.Exp3"), y = ("DBPct.Invivo"), 
          color = "Class", shape = 16, size = 2.5,alpha = 0.9, # Points color, shape and size
          palette = as.character(mycolors$Clr_list.1.18.),
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.method = "pearson",
          xlab = "Double bond %mol in primary culture #3", ylab = "Double bond %mol in vivo"
          )
a+ facet_wrap(~Age)+
  theme(text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))
ggsave(filename = "./Figure_Panels/Fig.1k.pdf", width = 6, height = 6, useDingbats = FALSE)
