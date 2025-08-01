
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ComplexHeatmap)

source("./Scripts/Function_scripts/Effect_size_functions.R")

##====Plotting function====
heatmap.plot <- function(plot.df, dataset){
  pdf(paste0("./Figure_Panels/EDFig.2.",dataset, ".pdf"), width =5, height =5)
  HM.2 <- Heatmap(plot.df, name = "z score", 
                  cluster_rows = TRUE,
                  column_title = dataset,
                  row_names_gp = gpar(fontsize = 3),
                  column_names_gp = gpar(fontsize = 5),
                  width = unit(5, "cm"),
                  height = unit(6.5, "cm"),
                  column_dend_height = unit(0.2, "cm"),
                  row_dend_width = unit(0.2, "cm"),
                  column_title_gp = gpar(fontsize = 5)
  )
  draw(HM.2)
  dev.off()
}

##====Primary culture #1====
load("./Output_Data/Exp1_Qui_sample_OvY_stats.Rdata")
load(file = "./Output_data/Spike-in_norm_Mednorm_all_373_lipids_back_to_raw_int.Rdata")

Exp1 <- raw_int %>% 
  rownames_to_column(var = "LipidIon")

Exp1.df <- ether.rename(Exp1)

q.sig <- Exp1.stat %>%
  rename("LipidIon" = "LipidIon...1") %>% 
  filter(., Padj<0.05) #95 lipids
sig.list <- ether.rename(q.sig)$LipidIon

Exp1.p <- Exp1.df %>% 
  filter(LipidIon %in% sig.list) %>% 
  select(matches("Quiescent_|LipidIon"))

z.sig.df <- Exp1.p %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  group_by(LipidIon) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(Conc))
  })

z.sig.p <- z.sig.df %>% 
  select(-Conc) %>% 
  pivot_wider(names_from = Sample, values_from = zscore) %>% 
  column_to_rownames(var = "LipidIon")

heatmap.plot(z.sig.p, "Primary_culture_#1")

##====Primary culture #2====
load(file = "./Output_Data/Exp2_Norm_Impt_backtoraw_conc612_lipids.Rdata")

Exp2.df <- raw_conc.exp2 %>% 
  rownames_to_column(var = "LipidIon") %>% 
  select(matches("_N|LipidIon"))

E2.n <- ether.rename(Exp2.df)

E2.l <- E2.n %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

E2.stat <- wilcox_stat(E2.l, Conc, LipidIon)

E2.n.sig <- E2.n %>% 
  filter(LipidIon %in% E2.stat$LipidIon...1[E2.stat$Wilcox_Pval<0.05])


z.sig.df2 <- E2.n.sig %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  group_by(LipidIon) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(Conc))
  })

z.sig.p2 <- z.sig.df2 %>% 
  select(-Conc) %>% 
  pivot_wider(names_from = Sample, values_from = zscore) %>% 
  column_to_rownames(var = "LipidIon")

heatmap.plot(z.sig.p2, "Primary_culture_#2")


##====Primary culture #3====
load(file = "./Output_data/Exp3_qNSC_quant.lipids.Rdata")

E3.n <- ether.rename(Exp3.conc.Q)

E3.l <- E3.n %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

E3.stat <- wilcox_stat(E3.l, Conc, LipidIon)

E3.n.sig <- E3.n %>% 
  filter(LipidIon %in% E3.stat$LipidIon...1[E3.stat$Wilcox_Pval<0.1])

z.sig.df3 <- E3.n.sig %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  group_by(LipidIon) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(Conc))
  })

z.sig.p3 <- z.sig.df3 %>% 
  select(-Conc) %>% 
  pivot_wider(names_from = Sample, values_from = zscore) %>% 
  column_to_rownames(var = "LipidIon")

heatmap.plot(z.sig.p3, "Also_Fig.1e_Primary_culture_#3")

##====In vivo====
load("./Output_Data/Invivo_Norm_Impt_log2_conc_121_lipids.Rdata")

Invivo <- 2^Invivo.Impt_norm_conc_all %>% 
  rownames_to_column(var = "LipidIon")

Invivo.n <- ether.rename(Invivo) 

Invivo.l <- Invivo.n %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

Invivo.stat <- wilcox_stat(Invivo.l, Conc, LipidIon)

Invivo.n.sig <- Invivo.n %>% 
  filter(LipidIon %in% Invivo.stat$LipidIon...1[Invivo.stat$Wilcox_Pval<0.1])

z.sig.df.invivo <- Invivo.n.sig %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  group_by(LipidIon) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(Conc))
  })

z.sig.p.invivo <- z.sig.df.invivo %>% 
  select(-Conc) %>% 
  pivot_wider(names_from = Sample, values_from = zscore) %>% 
  column_to_rownames(var = "LipidIon")

heatmap.plot(z.sig.p.invivo, "In_vivo")

##====Lipidyzer====
load("./Output_Data/Lipidyzer_qNSC_dup_TG_removed.Rdata")

ldz.l <- Ldz.Qui.nodup %>% 
  select(-ID_string) %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

ldz.stat <- wilcox_stat(ldz.l, Conc, LipidIon)

ldz.n.sig <- Ldz.Qui.nodup %>% 
  filter(LipidIon %in% ldz.stat$LipidIon...1[ldz.stat$Wilcox_Pval<0.1]) #65 lipids

z.sig.df.ldz <- ldz.n.sig %>%
  select(-ID_string) %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc") %>%
  group_by(LipidIon) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(Conc))
  })

z.sig.p.ldz <- z.sig.df.ldz %>% 
  select(-Conc) %>% 
  pivot_wider(names_from = Sample, values_from = zscore) %>% 
  mutate(., ID = str_split(LipidIon, "\\.")) %>% 
  mutate(Class = substr(LipidIon, 1, stri_locate_first(LipidIon, regex = "\\.")-1)) %>% 
  mutate(Class = ifelse(grepl("\\.O\\.|\\.P\\.", LipidIon), paste0(Class, "_Ether"), Class)) %>% 
  rowwise() %>% 
  mutate(LipidID = case_when(
    Class =="DAG" ~ paste0("DG(", ID[2], ":", ID[3], "_", ID[4], ":", ID[5], ")"), 
    Class =="PC" ~ paste0("PC(", ID[2], ":", ID[3], "_", ID[4], ":", ID[5], ")"),
    Class =="PE" ~ paste0("PE(", ID[2], ":", ID[3], "_", ID[4], ":", ID[5], ")"),
    Class =="PE_Ether" ~ paste0("PE(P-", ID[3], ":", ID[4], "_", ID[5], ":", ID[6], ")"),
    Class =="LPC" ~ paste0("LPC(", ID[2], ":", ID[3], ")"), 
    Class =="CER" ~ paste0("CER(", ID[2], ":", ID[3], ")"), 
    Class =="HCER" ~ paste0("HCER(", ID[2], ":", ID[3], ")"), 
    Class =="FFA" ~ paste0("FFA(", ID[2], ":", ID[3], ")"), 
    Class =="SM" ~ paste0("SM(", ID[2], ":", ID[3], ")"), 
    grepl("^TAG", LipidIon) ~ paste0("TG(", substr(ID[1], str_locate(ID[1], "G")+1, nchar(ID[1])), ":", ID[2], ")")
  )) %>% 
  ungroup() %>% 
  select(-c(Class, LipidIon, ID)) %>% 
  column_to_rownames(var = "LipidID")


heatmap.plot(z.sig.p.ldz, "Lipidyzer")
