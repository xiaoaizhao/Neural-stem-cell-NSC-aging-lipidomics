## Free fatty acid quantification in Primary Culture #2
rm(list=ls())
library(tidyverse)
library(ggpubr)

setwd(rstudioapi::getActiveProject())

## -------------------------------------------------------------------------------------------------------------------
FFA <- read.csv(file = "./Input_Data/FFA.csv", stringsAsFactors = F)

load(file = "./Output_Data/Exp2_Norm_factor_from_613lpds.Rdata") # load in normalization factor from complex lipid quantification
norm.df <- norm_mtx %>%
  rownames_to_column(., var = "Sample")

## -------------------------------------------------------------------------------------------------------------------
##normalize by median concentration of complex lipids first####
FFA.df <- FFA %>%
  rename(., "Sample" = "X") %>%
  filter(., !grepl("QC|PBS", Sample))

FFA.all <- left_join(FFA.df, norm.df, by = "Sample") %>%
  select(., -MedConc) %>%
  column_to_rownames(., var = "Sample")

FFA.norm <- sweep(FFA.all, 1, FFA.all$ind, "/")

FFA.norm.df <- FFA.norm %>%
  rownames_to_column(., var = "Sample")

##calculate concentration of all FFA based on Oleic Acid (d17)####
FFA.mtx <- FFA.norm.df %>%
  select(., -c(ind, D.17_Oleic_Acid)) %>%
  pivot_longer(-Sample, names_to = "Endo_FFA", values_to = "FFA_Int") %>%
  mutate(., Endo_FFA = str_replace(Endo_FFA, "\\.", ":")) %>%
  mutate(., Endo_FFA = str_replace(Endo_FFA, "X", ""))

IS_OA <- FFA.norm.df[, c("Sample", "D.17_Oleic_Acid")]

FFA.conc <- left_join(FFA.mtx, IS_OA, by = "Sample") %>%
  mutate(., FFA_conc = 0.5*FFA_Int/D.17_Oleic_Acid) %>%
  select(., c("Sample", "Endo_FFA", "FFA_conc")) %>%
  pivot_wider(., names_from = Endo_FFA, values_from = FFA_conc) 

save(FFA.conc, file = "./Output_Data/Exp2_FFA_concentration.rdata")



