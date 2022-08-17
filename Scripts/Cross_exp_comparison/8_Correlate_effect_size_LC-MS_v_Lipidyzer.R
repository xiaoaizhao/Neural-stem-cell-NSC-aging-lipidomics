
##Correlate effect sizes between lipidyzer and LC-MS data

##Steps:
#1. Load in effect size matrix of LC-MS and lipidyzer data
#2. Subset to only keep commonly detected lipids (115 lipids total)
#3. Make correlation plot between Age comparison and cell type comparison

##Results: very good correlation between lc-ms and lipidyzer both in age-related differences and cell type differences.
rm(list=ls())
library(tidyverse)
library(ggpubr)
setwd(rstudioapi::getActiveProject())

##Effect size between old vs. young qNSCs####=============================================================================
load("./Output_Data/Ef_Size_Lipid_Age_qNSC_InVitro_LC.Rdata")
load("./Output_Data/Ldz_Ef_Size_Age_qNSC.Rdata")
load("./Output_Data/LC_Lipidyzer_overlap_115_lipids_from_6_clas.Rdata")

LC_dyzer_overlap <- LC_dyzer_overlap %>%
  mutate(., LipidIon...2 = LipidIon.x) %>%
  mutate(., LipidIon...1 = LipidIon.y) 
InVitro.lpd.es.g <- InVitro.lpd.es.g %>%
  rename(., LipidIon...2 = LipidIon) %>%
  rename(., es_g_LC = es_g)

LC.filter <- inner_join(LC_dyzer_overlap, InVitro.lpd.es.g, by = "LipidIon...2") 

dyzer.effect.size <- Lipidyzer.Age.es.g %>%
  rename(., LipidIon...1 = Lipid) %>%
  rename(., es_g_Dyzer = es_g)
Age.corr.LC.Dyzer <- inner_join(LC.filter, dyzer.effect.size, by = "LipidIon...1")

ggscatter(Age.corr.LC.Dyzer, x = "es_g_Dyzer", y = "es_g_LC", 
          add = "reg.line", conf.int = TRUE, 
          alpha = 0.7,
          cor.coef = TRUE, 
          cor.method = "pearson",
          xlab = "Lipidyzer old vs. young qNSC", ylab = "LC-MS old vs. young qNSC",
          title = "Lipidyzer vs. LC-MS (effect size between age in qNSC)",
          font.label = c(12, "plain")
)
ggsave(filename = paste0("./Figure_Panels/Fig_S1h.pdf"), width = 5, height = 5,
       useDingbats=FALSE)

##Effect size between quiescent (qNSCs) vs. activated (aNSCs)####====================================================================================
load("./Output_Data/Ef_Size_Lipid_CellType_InVitro_LC.Rdata")
load("./Output_Data/Ldz_Ef_Size_celltype.Rdata")
load("./Output_Data/LC_Lipidyzer_overlap_115_lipids_from_6_clas.Rdata")

LC_dyzer_overlap <- LC_dyzer_overlap %>%
  mutate(., LipidIon...2 = LipidIon.x) %>%
  mutate(., LipidIon...1 = LipidIon.y) 

InVitro.cell.es.g <- Invitro.cell.es.g %>%
  rename(., LipidIon...2 = LipidIon) %>%
  rename(., es_g_LC = es_g)

LC.filter <- inner_join(LC_dyzer_overlap, InVitro.cell.es.g, by = "LipidIon...2") 

dyzer.effect.size <- Lipidyzer.cell.es.g %>%
  rename(., LipidIon...1 = Lipid) %>%
  rename(., es_g_Dyzer = es_g)
Cell.corr.LC.Dyzer <- inner_join(LC.filter, dyzer.effect.size, by = "LipidIon...1")

ggscatter(Cell.corr.LC.Dyzer, x = "es_g_Dyzer", y = "es_g_LC", 
          add = "reg.line", conf.int = TRUE, 
          alpha = 0.7,
          cor.coef = TRUE, 
          cor.method = "pearson",
          xlab = "Lipidyzer aNSC vs. qNSC", ylab = "LC-MS aNSC vs. qNSC",
          title = "Lipidyzer vs. LC-MS (effect size between cell type)",
          font.label = c(12, "plain")
)
ggsave(filename = paste0("./Figure_Panels/Fig_S1g.pdf"), width = 5, height = 5,
       useDingbats=FALSE)
