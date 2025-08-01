
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(stringi)
library(ggpubr)
source("./Scripts/Function_scripts/Effect_size_functions.R")

## ====Primary culture #3 data====
load("./Output_data/Exp3_backtoRAW_Norm_Impt_log2_conc_356_noconc_60_lipids.Rdata")

B3.LS.Cla<- raw_int.Exp3%>% 
  rownames_to_column(var = "LipidIon") 

Smp.Key <- read.csv("./Input_Data/Batch3_sample_key_forR.csv", stringsAsFactors = F)

Smp.Key.e <- Smp.Key %>% 
  filter(!grepl("PBS", Sample.Name))

Exp3.Q.lpd<- B3.LS.Cla %>% 
  rename_at(vars(matches(Smp.Key.e$Sample_ID)), ~Smp.Key.e$Sample.Name) %>% 
  filter(!grepl("Cholesterol\\+H\\-H2O_23_Positive|Ch\\+H\\-H2O", LipidIon)) %>% 
  select(-"O8_aNSC-A") %>% 
  select(matches("_qNSC-Q|LipidIon"))

E3.l <- ether.rename(Exp3.Q.lpd) %>% 
  pivot_longer(-LipidIon, values_to = "Conc", names_to = "Samples") %>% 
  mutate(Age = case_when(
    grepl("^Y", Samples) ~ "Young",
    grepl("^O", Samples) ~ "Old"
  )) %>% 
  group_by(LipidIon, Age) %>% 
  summarise(., MeanConc = mean(Conc))



## ====In vivo data====
load("./Output_Data/Invivo_Norm_Impt_log2_conc_121_lipids.Rdata")
Invivo.lpd <- 2^Invivo.Impt_norm_conc_all %>% 
  rownames_to_column(var = "LipidIon")

Invivo <- ether.rename(Invivo.lpd) 

Invivo.l <- Invivo %>% 
  ungroup() %>% 
  pivot_longer(-LipidIon, values_to = "Conc", names_to = "Samples") %>% 
  mutate(Age = case_when(
    grepl("^Y", Samples) ~ "Young",
    grepl("^O", Samples) ~ "Old"
  )) %>% 
  group_by(LipidIon, Age) %>% 
  summarise(., MeanConc = mean(Conc))


## ====Merge datasets====
Conc.E3invivo <- inner_join(E3.l, Invivo.l, by = c("Age", "LipidIon")) %>% 
  mutate(Class = ifelse(grepl("\\(",LipidIon), substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1),
                        LipidIon)) %>% 
  mutate(Log2Exp3 = log2(MeanConc.x)) %>% 
  mutate(Log2Invivo = log2(MeanConc.y))

Conc.E3invivo$Age <- factor(Conc.E3invivo$Age, levels = c("Young", "Old"))

mem <- c("PC", "PE", "PI", "PS", "Cholesterol")
TG <- c("TG")
Lyso <- c("LPE", "LPC")

df.mem <- Conc.E3invivo %>% 
  filter(Class %in% mem)

df.TG <- Conc.E3invivo %>% 
  filter(Class %in% TG)

df.Lyso <- Conc.E3invivo %>% 
  filter(Class %in% Lyso)

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")
mycolors <- lut.march25 %>% 
  filter(lipid_cat %in% df.mem$Class)

mycolors$lipid_cat <- factor(mycolors$lipid_cat)

mycolors <- mycolors %>% 
  arrange(lipid_cat)

## ====Plot phospholipids and cholesterol====
df.mem$Class <- factor(df.mem$Class, levels = levels(mycolors$lipid_cat))

a <- ggscatter(df.mem, x = ("Log2Exp3"), y = ("Log2Invivo"), 
          color = "Class", shape = 16, size = 2.5,alpha = 0.9, # Points color, shape and size
          palette = as.character(mycolors$Clr_list.1.18.),
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.method = "pearson",
          xlab = "Log2 lipid concentration (uM) in primary culture #3", ylab = "Log2 lipid concentration (uM) in vivo",
          title = "Phospholipids and cholesterol")
a+ facet_wrap(~Age)+
  theme(text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))
ggsave(filename = "./Figure_Panels/Fig.1i.pdf", width = 6, height = 6, useDingbats = FALSE)


mycolors <- lut.march25 %>% 
  filter(lipid_cat %in% df.TG$Class)

mycolors$lipid_cat <- factor(mycolors$lipid_cat)

mycolors <- mycolors %>% 
  arrange(lipid_cat)

## ====Plot glycerol lipids====
df.TG$Class <- factor(df.TG$Class, levels = levels(mycolors$lipid_cat))

a <- ggscatter(df.TG, x = ("Log2Exp3"), y = ("Log2Invivo"), 
          color = "Class", shape = 16, size = 2.5,alpha = 0.9, # Points color, shape and size
          palette = as.character(mycolors$Clr_list.1.18.),
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.method = "pearson",
          xlab = "Log2 lipid concentration (uM) in primary culture #3", ylab = "Log2 lipid concentration (uM) in vivo",
          title = "Glycerol lipids")
a+ facet_wrap(~Age)+
  theme(text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))
ggsave(filename = "./Figure_Panels/EDFig.5d.pdf", width = 6, height = 6, useDingbats = FALSE)

## ====Plot lysophospholipids ====
mycolors <- lut.march25 %>% 
  filter(lipid_cat %in% df.Lyso$Class)

mycolors$lipid_cat <- factor(mycolors$lipid_cat)

mycolors <- mycolors %>% 
  arrange(lipid_cat)

df.Lyso$Class <- factor(df.Lyso$Class, levels = levels(mycolors$lipid_cat))

a <- ggscatter(df.Lyso, x = ("Log2Exp3"), y = ("Log2Invivo"), 
          color = "Class", shape = 16, size = 2.5,alpha = 0.9, # Points color, shape and size
          palette = as.character(mycolors$Clr_list.1.18.),
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.method = "pearson",
          xlab = "Log2 lipid concentration (uM) in primary culture #3", ylab = "Log2 lipid concentration (uM) in vivo",
          title = "Lysophospholipids")
a+ facet_wrap(~Age)+
  theme(text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))
ggsave(filename = "./Figure_Panels/EDFig.5c.pdf", width = 6, height = 6, useDingbats = FALSE)
