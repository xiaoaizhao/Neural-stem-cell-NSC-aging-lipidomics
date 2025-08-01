

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
source("./Scripts/Function_scripts/Effect_size_functions.R")

load("./Output_Data/Lipidyzer_aNSC_for.LC-MS.ovlp.Rdata")
load("./Output_Data/Lipidyzer_qNSC_for.LC-MS.ovlp.Rdata")
load("./Output_Data/Exp3_lipid.conc.reformat.for.lipidyzer.overlap.Rdata")

## ==== Primary culture #3, get mean lipid concentration in each age group====
E3.mean <-  E3.lpd.conc.fmt %>% 
  select(-LipidIon1) %>% 
  pivot_longer(-c(LipidIon, ID_string), values_to = "Concentration", names_to = "Samples") %>% 
  mutate(Age = ifelse(grepl("^Y", Samples), "Young", "Old")) %>% 
  mutate(CellType = ifelse(grepl("_qNSC-Q", Samples), "qNSC", "aNSC")) %>% 
  group_by(LipidIon, ID_string, Age, CellType) %>% 
  summarise(log2.Exp3.mean = log2(mean(Concentration)))

E3.Q.mean <- E3.mean %>% 
  filter(CellType == "qNSC")

E3.A.mean <- E3.mean %>% 
  filter(CellType == "aNSC")

## ==== Lipidyzer qNSC, get mean lipid concentration in each age group====
Ldz.Q.mean <-  Ldz.Qui.nodup %>% 
  pivot_longer(-c(LipidIon, ID_string), values_to = "Concentration", names_to = "Samples") %>% 
  mutate(Age = ifelse(grepl("^Y", Samples), "Young", "Old")) %>% 
  group_by(LipidIon, ID_string, Age) %>% 
  summarise(log2.Ldz.Q.mean = log2(mean(Concentration)))

## ==== Overlap in qNSCs ====
qNSC.ovlp <- inner_join(E3.Q.mean, Ldz.Q.mean, by = c("ID_string", "Age")) %>% 
  mutate(Class = substr(LipidIon.x, 1, str_locate(LipidIon.x, "\\(")-1))

qNSC.ovlp$Age <- factor(qNSC.ovlp$Age, levels = c("Young", "Old"))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")
mycolors <- lut.march25 %>% 
  filter(lipid_cat %in% qNSC.ovlp$Class)

mycolors$lipid_cat <- factor(mycolors$lipid_cat)

mycolors <- mycolors %>% 
  arrange(lipid_cat)

qNSC.ovlp$Class <- factor(qNSC.ovlp$Class, levels = levels(mycolors$lipid_cat))
qNSC.ovlp$Age <- factor(qNSC.ovlp$Age, levels = c("Young", "Old"))

a <- ggscatter(qNSC.ovlp, x = ("log2.Exp3.mean"), y = ("log2.Ldz.Q.mean"), 
               color = "Class", shape = 16, size = 2.5,alpha = 0.9, # Points color, shape and size
               palette = as.character(mycolors$Clr_list.1.18.),
               add = "reg.line",  # Add regressin line
               add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
               conf.int = TRUE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
               cor.method = "pearson",
               xlab = "Log2 uM concentration (Primary culture #3)", ylab = "Log2 nmol/g concentration (Lipidyzer)",
               title = "qNSC - Primary culture #3 vs. Lipidzyer")
a+ facet_wrap(~Age)+
  theme(text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))
ggsave(filename = "./Figure_Panels/EDFig.3d.pdf", width = 6, height = 6, useDingbats = FALSE)


## ==== Lipidyzer aNSC, get mean lipid concentration in each age group====
Ldz.A.mean <-  Ldz.Act.nodup %>% 
  pivot_longer(-c(LipidIon, ID_string), values_to = "Concentration", names_to = "Samples") %>% 
  mutate(Age = ifelse(grepl("^Y", Samples), "Young", "Old")) %>% 
  group_by(LipidIon, ID_string, Age) %>% 
  summarise(log2.Ldz.A.mean = log2(mean(Concentration)))

## ==== Overlap in aNSCs ====
aNSC.ovlp <- inner_join(E3.A.mean, Ldz.A.mean, by = c("ID_string", "Age")) %>% 
  mutate(Class = substr(LipidIon.x, 1, str_locate(LipidIon.x, "\\(")-1))

aNSC.ovlp$Age <- factor(aNSC.ovlp$Age, levels = c("Young", "Old"))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")
mycolors <- lut.march25 %>% 
  filter(lipid_cat %in% aNSC.ovlp$Class)

mycolors$lipid_cat <- factor(mycolors$lipid_cat)

mycolors <- mycolors %>% 
  arrange(lipid_cat)

aNSC.ovlp$Class <- factor(aNSC.ovlp$Class, levels = levels(mycolors$lipid_cat))
aNSC.ovlp$Age <- factor(aNSC.ovlp$Age, levels = c("Young", "Old"))

a <- ggscatter(aNSC.ovlp, x = ("log2.Exp3.mean"), y = ("log2.Ldz.A.mean"), 
               color = "Class", shape = 16, size = 2.5,alpha = 0.9, # Points color, shape and size
               palette = as.character(mycolors$Clr_list.1.18.),
               add = "reg.line",  # Add regressin line
               add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
               conf.int = TRUE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
               cor.method = "pearson",
               xlab = "Log2 uM concentration (Primary culture #3)", ylab = "Log2 nmol/g concentration (Lipidyzer)",
               title = "aNSC - Primary culture #3 vs. Lipidzyer")
a+ facet_wrap(~Age)+
  theme(text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))
ggsave(filename = "./Figure_Panels/EDFig.3e.pdf", width = 6, height = 6, useDingbats = FALSE)