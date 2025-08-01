
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(patchwork)
library(stringi)

source("./Scripts/Function_scripts/Effect_size_functions.R")

DB.sig.fraction <- function(df, dataset.name) {
load(file = "./Output_Data/DB.2outof3LC.invitro.Sig.features.Rdata")
DB.sig.ls <- unique(p.c2.in.3lc.df.sig$Cla_DB)

df1 <- df %>% 
  mutate(Cla_DB = paste0(Class, DB_num)) %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old")) %>% 
  group_by(Cla_DB, Age) %>% 
  summarise(mean.DB.conc = mean(Sum_DB)) %>% 
  filter(Cla_DB %in% DB.sig.ls)

df2 <- df1 %>%
  group_by(Age) %>% 
  summarise(Sum.DB.conc = sum(mean.DB.conc))

df3 <- left_join(df1, df2, by = "Age") %>% 
  mutate(Pct_conc = (mean.DB.conc / Sum.DB.conc) * 100) %>% 
  mutate(Class = substr(Cla_DB, 1, str_locate(Cla_DB, ":")-1))


load("./Output_Data/Class_col_list_paper.order_031725.Rdata")

mycolors <- lut.march25 %>% 
  filter(lipid_cat %in%df3$Class)

df3$Class <- factor(df3$Class, levels = c(mycolors$lipid_cat))

df3 <- df3 %>% 
    arrange(Class)

df3$Age <- factor(df3$Age, levels = c("Young", "Old"))

b <- ggplot(df3, aes(x=Age, y=Pct_conc))+
  geom_bar(aes(fill=Class),
  position="stack", stat="identity")+
    labs(title  = paste0("DB sig ftr proportion", dataset.name))+
    ylab("Proportional concentration (mol%) ")+
    xlab("")+
    scale_fill_manual(values = as.character(mycolors$Clr_list.1.18.)) +
    theme(legend.position= "right")+
    guides(fill=guide_legend(title="Class")) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"))

b

ggsave(filename = paste0("./Figure_Panels/EDFig.4d.",dataset.name, ".pdf"), width = 5, height = 5, useDingbats = FALSE)
}  


##====Primary culture #1====
load("./Output_Data/Qui_NSC_DB_PCT_by_Class_2017LC-MS.Rdata")
DB.sig.fraction(Qui_LC2017_df,  "Primary_culture_#1")

##====Primary culture #2====
load("./Output_Data/Exp2_CONC.DB_PCT_all_samples.Rdata")
Ctrl.exp2 <- Exp2_CONC.DB %>% 
  filter(KO == "_N")
DB.sig.fraction(Ctrl.exp2,  "Primary_culture_#2")

##====Primary culture #3====
load("./Output_Data/Exp3.Qui_CONC.DB_PCT_by_Class.Rdata")
DB.sig.fraction(Qui_CONC.DB.E3,  "Primary_culture_#3")
