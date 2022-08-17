
## Effect size calculation using Hedge's g on primary culture #2 data
## Effect size of individual lipids
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")

#### Effect size on individual lipid species #########################################################################################
##calculate effect size on primary culture #2 data####

load(file = "./Output_Data/Exp2_Norm_Impt_backtoraw_all693_lipids.Rdata") ##load pre-processed matrix with lipid concentration
all.lipid <- raw_int.exp2 %>%  
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") 

all.lipid.df <- all.lipid %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old")) %>%
  mutate(., KO = substr(Sample, nchar(Sample), nchar(Sample)))

##effect size on age difference at individual lipid level within each KO condition
KO.list <- unique(all.lipid.df$KO)

EfSize.df = list()
for (KO.name in KO.list) {
    KO.df<- all.lipid.df %>%
      filter(., KO == KO.name)
    EfSize.df[[KO.name]] <- es.g.func(KO.df, LipidIon, Age, Conc_Int, Sample) %>%
      mutate(., LipidIon = str_replace_all(LipidIon, "/", "_")) %>%
      mutate(., KO = KO.name)
}

Exp2.lpd.es.g.allKO <- bind_rows(EfSize.df)
save(Exp2.lpd.es.g.allKO, file = "./Output_Data/Ef_Size_Lipid_Exp2_all_KO.Rdata")

#### Effect size on double bond composition #########################################################################################
#load double bond %mol calculated from previous scripts
load("./Output_Data/Exp2_DB_PCT_all_samples.Rdata")
Exp2_DB <- Exp2_DB %>%
  mutate(., Cla_DB = paste0(Class, DB_num)) %>%
  mutate(., KO = substr(Sample, nchar(Sample), nchar(Sample)))

##effect size on age difference at double bond composition level within each KO condition
KO.list <- unique(Exp2_DB$KO)

EfSize.df.DB = list()
for (KO.name in KO.list) {
  KO.df<- Exp2_DB %>%
    filter(., KO == KO.name)
  EfSize.df.DB[[KO.name]] <- es.g.func(KO.df, Cla_DB, Age, DB_Pct, Sample) %>%
    mutate(., KO = KO.name)
}

Exp2.DB.es.g.allKO <- bind_rows(EfSize.df.DB)
save(Exp2.DB.es.g.allKO, file = "./Output_Data/Ef_Size_DB_PCT_Exp2_all_KO.Rdata")

