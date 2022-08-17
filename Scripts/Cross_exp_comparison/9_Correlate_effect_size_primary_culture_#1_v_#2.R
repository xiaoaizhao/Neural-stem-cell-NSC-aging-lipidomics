
##Correlation analysis between primary culture (2017 LC-MS) data and primary culture #2 (KO control data in 2019)
##Using Hedge's g effect size
##Subset to only include control sample from primary culture #2 (non-targeting control) and then generate correlation plot.
rm(list=ls())
library(tidyverse)
library(ggpubr)
setwd(rstudioapi::getActiveProject())

load("./Output_Data/Ef_Size_Lipid_Age_qNSC_InVitro_LC.Rdata")
load("./Output_Data/Ef_Size_Lipid_Exp2_all_KO.Rdata")

Set2.NTG <- Exp2.lpd.es.g.allKO %>%
  filter(., KO == "N")
Eff_size_all <- inner_join(Set2.NTG, InVitro.lpd.es.g, by = "LipidIon") %>%
  select(., matches("LipidIon|es.g")) %>%
  rename(., "Effect_size_KO"= "es_g.x") %>%
  rename(., "Effect_size_invitro"= "es_g.y") 

####correlation between two datasets####

ggscatter(Eff_size_all, x = "Effect_size_KO", y = "Effect_size_invitro", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, 
          cor.method = "pearson",
          xlab = "Primary NSC Culture #2", ylab = "Primary NSC Culture",
          title = "Effect Size Correlation",
          alpha = 0.7,
          font.label = c(12, "plain"),
)
ggsave(filename = paste0("./Figure_Panels/Fig_S1d.pdf"), height = 5, width = 5,
       useDingbats=FALSE)
