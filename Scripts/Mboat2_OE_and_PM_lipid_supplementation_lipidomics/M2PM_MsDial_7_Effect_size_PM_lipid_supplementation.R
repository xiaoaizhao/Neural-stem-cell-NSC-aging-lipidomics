### Effect size between control and GPMV-supplemented samples
#(pos ES -> higher in Treatment)
#(neg ES -> higher in Control)
rm(list=ls())
library(tidyverse)
source("./Scripts/Function_scripts/Effect_size_functions.R")

# === Effect size calculation on individual lipids====
load(file = "./Output_data/PM_supp_LIPID.Rdata")

## Remove outlier sample: XZ_45 - Y5_YlpdP5
PM <- PM.sup.Lipid %>% 
  select(-Y5_YlpdP5) %>% 
  rownames_to_column(var = "LipidIon") %>% 
  pivot_longer(-LipidIon, names_to = "Samples", values_to = "Conc_Int") %>% 
  mutate(., Age = case_when(
    grepl("^Y", Samples) ~ "Young",
    grepl("^O", Samples) ~ "Old"
  )) %>% 
  mutate(., Condition = case_when(
    grepl("_meth", Samples) ~ "Control",
    grepl("_Ylpd", Samples) ~ "Treatment"
  ))


EfSize.df = list()
Ages <- c("Young", "Old")
for (age in Ages) {
  t.df<- PM %>%
    filter(., Age == age)
  EfSize.df[[age]] <- es.g.treat.func(t.df, LipidIon, Condition, Conc_Int, Samples) %>%
    mutate(., LipidIon = str_replace_all(LipidIon, "/", "_")) %>%
    mutate(., Age = age)
}
PM.supp.ef.sz<- bind_rows(EfSize.df)

save(PM.supp.ef.sz, file = "./Output_Data/PM_lipid_supplementation_by_treatment_all_ages.Rdata")

