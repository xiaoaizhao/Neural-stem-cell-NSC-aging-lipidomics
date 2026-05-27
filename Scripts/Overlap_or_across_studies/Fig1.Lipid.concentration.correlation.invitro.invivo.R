

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(stringi)
library(ggpubr)
library(patchwork)
source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")

## ====In vitro====
load("./Output_data/Invitro_conc.lipid.qNSC.Rdata")

Invitro.l <- MsD.lpd.rmv.abc(Invitro.q.conc) %>% 
  pivot_longer(-LipidIon, values_to = "Conc", names_to = "Samples") %>% 
  mutate(Age = case_when(
    grepl("^Y", Samples) ~ "Young",
    grepl("^O", Samples) ~ "Old"
  )) %>% 
  group_by(LipidIon, Age) %>% 
  summarise(., MeanConc = mean(Conc))


## ====In vivo data====
load("./Output_Data/Invivo.MsD.29.lipid.raw.conc.Rdata")
Invivo.df <- Invivo.raw.conc.29.good.peak 
Invivo.l <- MsD.lpd.rmv.abc(Invivo.df) %>% 
  ungroup() %>% 
  pivot_longer(-LipidIon, values_to = "Conc", names_to = "Samples") %>% 
  mutate(Age = case_when(
    grepl("^Y", Samples) ~ "Young",
    grepl("^O", Samples) ~ "Old"
  )) %>% 
  group_by(LipidIon, Age) %>% 
  summarise(., MeanConc = mean(Conc))


## ====Merge datasets====
Conc.Invitro.invivo <- inner_join(Invitro.l, Invivo.l, by = c("Age", "LipidIon")) %>% 
  mutate(Class = ifelse(grepl("\\(",LipidIon), substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1),
                        LipidIon)) %>% 
  mutate(Log2Invitro = log2(MeanConc.x)) %>% 
  mutate(Log2Invivo = log2(MeanConc.y))

Conc.Invitro.invivo$Age <- factor(Conc.Invitro.invivo$Age, levels = c("Young", "Old"))

mem <- c("PC", "PE", "PI", "PS", "Cholesterol", "Cer")
Lyso_TG <- c("LPE", "LPC", "TG")

df.mem <- Conc.Invitro.invivo %>%
  filter(Class %in% mem)

df.Lyso_TG <- Conc.Invitro.invivo %>%
  filter(Class %in% Lyso_TG)

# ## ==== Plot  ====
load("./Output_Data/Class_col_list_paper.order_031725.Rdata")


## ====Plot lysophospholipids and glycerol lipids====
mycolors <- lut.march25 %>% 
  filter(lipid_cat %in% df.Lyso_TG$Class)

mycolors$lipid_cat <- factor(mycolors$lipid_cat)

mycolors <- mycolors %>% 
  arrange(lipid_cat)


df.Lyso_TG$Class <- factor(df.Lyso_TG$Class, levels = levels(mycolors$lipid_cat))


a <- ggscatter(df.Lyso_TG, x = ("Log2Invitro"), y = ("Log2Invivo"), 
          color = "Class", shape = 16, size = 2.5,alpha = 0.9, # Points color, shape and size
          palette = as.character(mycolors$Clr_list.1.18.),
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.method = "pearson",
          xlab = "Log2 lipid concentration (uM) In vitro", ylab = "Log2 lipid concentration (uM) In vivo",
          title = "Lysophospholipids and glycerol lipids") + 
  facet_wrap(~Age)+
  theme(text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))
a
ggsave(filename = "./Figure_Panels/fig.S6D.pdf", width = 6, height = 6, useDingbats = FALSE)

## ====Plot phospholipids and cholesterol====
mycolors <- lut.march25 %>% 
  filter(lipid_cat %in% df.mem$Class)

mycolors$lipid_cat <- factor(mycolors$lipid_cat)

mycolors <- mycolors %>% 
  arrange(lipid_cat)


df.mem$Class <- factor(df.mem$Class, levels = levels(mycolors$lipid_cat))


b <- ggscatter(df.mem, x = ("Log2Invitro"), y = ("Log2Invivo"), 
               color = "Class", shape = 16, size = 2.5,alpha = 0.9, # Points color, shape and size
               palette = as.character(mycolors$Clr_list.1.18.),
               add = "reg.line",  # Add regressin line
               add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
               conf.int = TRUE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
               cor.method = "pearson",
               xlab = "Log2 lipid concentration (uM) In vitro", ylab = "Log2 lipid concentration (uM) In vivo",
               title = "Phospho- and sphingolipids") + 
  facet_wrap(~Age)+
  theme(text = element_text(size = 12, face = "plain", colour = "black"))+
  theme(axis.text = element_text(size = 12, face = "plain", colour = "black"))

b
ggsave(filename = "./Figure_Panels/fig.S6C.pdf", width = 6, height = 6, useDingbats = FALSE)
