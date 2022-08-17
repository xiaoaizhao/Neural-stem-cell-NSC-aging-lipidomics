##Class fold change comparison between whole cell (Primary culture #1) vs. GPMV
##Steps:
#1. Get a list of common lipid classes detected in both whole cell and GPMV dataset 
#2. Remove intracellular lipid classes from lipid class list
#3. Calculate log2 fold change on sum concentration of each class between old and young sample
#4. Plot log2 fold change between datasets on common lipid classes exclude intracellular lipid classes.
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(ggrepel)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")


load("./Output_Data/GPMV_ClassSum_565_lipids.Rdata")
load("./Output_Data/Qui_NSC_ClassSum_Int_2017LC-MS.Rdata")

#### Common lipid classes####
cmn.ls <- intersect(classsum$Class, GPMV.classsum$Class) #17 common lipid classes
##Intracellular classes
intra.ls <- c("AcCa", "CL", "Co", "ZyE", "TG", "DG")

mem.ls <- cmn.ls[!cmn.ls %in% intra.ls]

#### Log2 fold change between young and old in whole cell data (Primary Culture #1)####
WC.df <- classsum %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old")) %>%
  group_by(Age)
WC.stat <- wilcox_stat(WC.df, Class_sum, Class) %>%
  filter(., Class...1 %in% mem.ls) %>%
  rename_all(paste0, "_Whole_Cell")

#### Log2 fold change between young and old in GPMV data####
GPMV.df <- GPMV.classsum %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old")) %>%
  group_by(Age)
GPMV.stat <- wilcox_stat(GPMV.df, Class_sum, Class) %>%
  filter(., Class...1 %in% mem.ls) %>%
  rename_all(paste0, "_GPMV")

FC.cmb <- bind_cols(WC.stat, GPMV.stat)

#### Plot ####
ggscatter(FC.cmb, x = "log2FC_OvY_GPMV", y = "log2FC_OvY_Whole_Cell", alpha = 0.8,
          cor.method = "pearson",
          xlab = "GPMV (Log2 Old/Young)", ylab = "Whole Cell (Log2 Old/Young)",
          title = "Class sum change with age (GPMV vs. Whole Cell) ",
          label = "Class...1_Whole_Cell", repel = TRUE
)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")
ggsave(filename = paste0("./Figure_Panels/Fig_S5c.pdf"), height = 5, width = 5, useDingbats=FALSE)
