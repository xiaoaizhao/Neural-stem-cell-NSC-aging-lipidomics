
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(stringi)
library(ggpubr)
source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_data/Exp3_qNSC_quant.lipids.Rdata")
Exp3.Q <- ether.rename(Exp3.conc.Q) 
E3.l <- Exp3.Q %>% 
  pivot_longer(-LipidIon, values_to = "Conc", names_to = "Samples") %>% 
  mutate(Age = case_when(
    grepl("^Y", Samples) ~ "Young",
    grepl("^O", Samples) ~ "Old"
  )) %>% 
  group_by(LipidIon, Age) %>% 
  summarise(., MeanConc = mean(Conc))

load("./Output_Data/Exp2_Norm_Impt_backtoraw_conc612_lipids.Rdata")
Exp2.lpd <- raw_conc.exp2 %>% 
  rownames_to_column(var = "LipidIon")

Exp2.l <- ether.rename(Exp2.lpd) %>% 
  pivot_longer(-LipidIon, values_to = "Conc", names_to = "Samples") %>% 
  mutate(., KO = substr(Samples, nchar(Samples)-1, nchar(Samples))) %>% 
  filter(KO == "_N") %>% 
  mutate(Age = case_when(
    grepl("^Y", Samples) ~ "Young",
    grepl("^O", Samples) ~ "Old"
  )) %>% 
  group_by(LipidIon, Age) %>% 
  summarise(., MeanConc = mean(Conc))


all.cor <- inner_join(Exp2.l, E3.l, by = c("LipidIon", "Age"))%>% 
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>% 
  mutate(Class = ifelse(grepl("^Cholesterol", LipidIon), "Cholesterol", Class)) %>% 
  mutate(., log2Exp3 = log2(MeanConc.y)) %>% 
  mutate(., log2Exp2 = log2(MeanConc.x))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")
mycolors <- lut.march25 %>%
  filter(lipid_cat %in% all.cor$Class)

mycolors$lipid_cat <- factor(mycolors$lipid_cat)

mycolors <- mycolors %>%
  arrange(lipid_cat)

all.cor$Class <- factor(all.cor$Class, levels = levels(mycolors$lipid_cat))
all.cor$Age <- factor(all.cor$Age, levels = c("Young", "Old"))

a <- ggscatter(all.cor, x = ("log2Exp3"), y = ("log2Exp2"), 
          color = "Class", shape = 16, size = 2.5,alpha = 0.9, # Points color, shape and size
          palette = as.character(mycolors$Clr_list.1.18.),
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.method = "pearson",
          xlab = "Log2 uM concentration (Primary culture #3)", ylab = "Log2 uM concentration (Primary culture #2)"
         )
a + facet_wrap(~Age)
ggsave(filename = "./Figure_Panels/EDFig.1d.pdf", width = 6, height = 6, useDingbats = FALSE)
