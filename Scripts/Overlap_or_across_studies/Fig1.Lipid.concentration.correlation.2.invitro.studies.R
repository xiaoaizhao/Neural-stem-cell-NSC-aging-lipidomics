### ----------------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(ggthemes)
library("scales")
library(ggpubr)
library(ggbeeswarm)

## -------------------------------------------------------------------------------------------------------------------------
rm(list = ls())
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Invitro_MsD.Norm_Impt_log2_conc_404_lipids.Rdata")
load("./Output_Data/Exp2_MsD.Norm_Impt_log2_conc_278_lipids.Rdata")

Invitro.df <- Invitro.log2.impt.norm.conc.MsD %>% 
  rownames_to_column(var = "LipidIon") 

Invitro <- MsD.lpd.rmv.abc(Invitro.df) %>% 
  select(c(LipidIon, contains("_qNSC-Q"))) %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

Invitro.age.mean <- Invitro %>% 
  group_by(LipidIon, Age) %>% 
  summarise(Conc_avg = mean(Conc)) %>% 
  rename("Log2_Invitro_Conc_avg" = "Conc_avg")

E2.df <- Exp2.log2.impt.norm.conc.MsD %>% 
  rownames_to_column(var = "LipidIon") %>% 
  select(LipidIon, contains("_N"))

E2 <- MsD.lpd.rmv.abc(E2.df) %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

E2.age.mean <- E2 %>% 
  group_by(LipidIon, Age) %>% 
  summarise(Conc_avg = mean(Conc)) %>% 
  rename("Log2_Invitro.E2_Conc_avg" = "Conc_avg")


## -------------------------------------------------------------------------------------------------------------------------
Conc.corr.2.invitro <- inner_join(Invitro.age.mean, E2.age.mean, by = c("LipidIon", "Age")) #350

Conc.corr.2.invitro$Age <- factor(Conc.corr.2.invitro$Age, levels = c("Young", "Old"))


## -------------------------------------------------------------------------------------------------------------------------
Conc.corr.2.invitro <- Conc.corr.2.invitro %>% 
  mutate(Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")
mycolors <- lut.march25 %>%
  filter(lipid_cat %in% Conc.corr.2.invitro$Class)

mycolors$lipid_cat <- factor(mycolors$lipid_cat)

mycolors <- mycolors %>%
  arrange(lipid_cat)

Conc.corr.2.invitro$Class <- factor(Conc.corr.2.invitro$Class, levels = levels(mycolors$lipid_cat))
Conc.corr.2.invitro$Age <- factor(Conc.corr.2.invitro$Age, levels = c("Young", "Old"))

a <- ggscatter(Conc.corr.2.invitro, x = ("Log2_Invitro_Conc_avg"), y = ("Log2_Invitro.E2_Conc_avg"), 
          color = "Class", shape = 16, size = 2.5,alpha = 0.9, # Points color, shape and size
          palette = as.character(mycolors$Clr_list.1.18.),
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.method = "pearson",
          xlab = "Log2 uM concentration (In vitro)", ylab = "Log2 uM concentration (In vitro Experiment #2)"
         )
a + facet_wrap(~Age)
ggsave(filename = "./Figure_Panels/EDFig.1d.pdf", width = 6, height = 6, useDingbats = FALSE)

