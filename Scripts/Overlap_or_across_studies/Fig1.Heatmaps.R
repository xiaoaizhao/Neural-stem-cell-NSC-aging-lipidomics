### Heatmaps
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ComplexHeatmap)

source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")
##====Plotting function====
heatmap.plot <- function(plot.df, dataset){
  pdf(paste0("./Figure_Panels/fig.S2.",dataset, ".pdf"), width =5, height =5)
  set.seed(12345)
  HM.2 <- Heatmap(plot.df, name = "z score", 
                  cluster_rows = TRUE,
                  column_title = dataset,
                  row_names_gp = gpar(fontsize = 3),
                  column_names_gp = gpar(fontsize = 5),
                  width = unit(5, "cm"),
                  height = unit(6.5, "cm"),
                  column_dend_height = unit(0.2, "cm"),
                  row_dend_width = unit(0.2, "cm"),
                  column_title_gp = gpar(fontsize = 5),
                  clustering_distance_rows = "euclidean",
                  clustering_method_rows = "average"
  )
  draw(HM.2)
  dev.off()
}

##====In vitro====
load(file = "./Output_data/Invitro_qNSCs_440_conc_lipid.MSD.format.Rdata")

Invitro.l <- E3MsD.frmt.q.conc %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

Invitro.stat <- wilcox_stat(Invitro.l, Conc, LipidIon)

Invitro.n.sig <- E3MsD.frmt.q.conc %>% 
  filter(LipidIon %in% Invitro.stat$LipidIon...1[Invitro.stat$Wilcox_Pval<0.1]) 

Invitro.n.sig.n <- MsD.lpd.rmv.abc(Invitro.n.sig) 

check.invitro.df <- Invitro.n.sig.n %>% 
  group_by(LipidIon) %>% 
  tally()#each lipid only appear once

z.sig.df3 <- Invitro.n.sig.n %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  group_by(LipidIon) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(Conc))
  }) %>% 
  ungroup() 


z.sig.p3 <- z.sig.df3 %>% 
  select(-Conc) %>% 
  pivot_wider(names_from = Sample, values_from = zscore) %>% 
  column_to_rownames(var = "LipidIon")

heatmap.plot(z.sig.p3, "Also_Fig.1d_Invitro")

##====In vitro Experiment #2====
load(file = "./Output_Data/Exp2_278_conc_lipid.MSD.format.Rdata")

E2.n <- MsD.lpd.format(E2MsD.frmt) %>% 
  select(matches("_N|LipidIon"))

E2.l <- E2.n %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old")) %>% 
  filter(!grepl("PE\\(P-17:0_18:1\\)b", LipidIon)) 

E2.stat <- wilcox_stat(E2.l, Conc, LipidIon)

E2.n.sig <- E2.n %>% 
  filter(LipidIon %in% E2.stat$LipidIon...1[E2.stat$Padj<0.1]) 

E2.n.sig.n <- MsD.lpd.rmv.abc(E2.n.sig) 

check.2.df <- E2.n.sig.n %>% 
  group_by(LipidIon) %>% 
  tally() #each lipid only appear once

z.sig.df2 <- E2.n.sig.n %>%
  group_by(LipidIon) %>% 
  group_modify(~{
    .x %>% 
      slice_sample(n = 1)
  }) %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  group_by(LipidIon) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(Conc))
  })

z.sig.p2 <- z.sig.df2%>% 
  select(-Conc) %>% 
  pivot_wider(names_from = Sample, values_from = zscore) %>% 
  column_to_rownames(var = "LipidIon")  

heatmap.plot(z.sig.p2, "In vitro Experiment #2")


##====In vivo====
load("./Output_Data/Invivo_MsD.Norm_Impt_log2_conc_29.lipids.Rdata")
Invivo.df <- Invivo.log2.impt.norm.conc.MsD.good.peak %>% 
  rownames_to_column(var = "LipidIon")
Invivo.n <- MsD.lpd.format(Invivo.df) 

Invivo.l <- Invivo.n %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

z.sig.df.invivo <- Invivo.l %>%
  ungroup() %>% 
  group_by(LipidIon) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(Conc))
  })

z.sig.p.invivo <- z.sig.df.invivo %>% 
  select(-c(Conc, Age)) %>% 
  pivot_wider(names_from = Sample, values_from = zscore) %>% 
  column_to_rownames(var = "LipidIon")

heatmap.plot(z.sig.p.invivo, "In_vivo")

