### Processing script 6, Primary culture #3
#1. Calculate concentration (remove IS labelled from endogenous data)
#2. Adjust concentration and intensity for O8-aNSC sample
setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(reshape2)
library(stringi)

load("./Output_data/IS_Exp3_w_conc.Rdata")
load("./Output_data/Exp3_LipMol_Ion+FA_clean_414.Rdata")
load("./Output_Data/IS_Exp_3.Rdata")

# Add two endogenouse cholesterol peaks from Xcalibur to the data matrix
## Organize dataframe from Xcalibur
Xcal.chol <- IS.batch3 %>% 
  select(., -c(matches("d7|d9")))

Xcal.chol.df <- as.data.frame(t(Xcal.chol))

colnames(Xcal.chol.df) <- Xcal.chol.df[1,]

Xcal.chol.df <- Xcal.chol.df[-1, ] %>% 
  select(-`9_1_Solvent_Blank_1`) %>% 
  rownames_to_column(., var = "LipidIon") %>% 
  mutate(., Class = "Cholesterol") %>% 
  mutate(., Ion.m = "+H-H2O") %>% 
  mutate_at(vars(matches("QC|XZ")), as.numeric)
 
## Append to data from LipidSearch
Exp3.all <- Exp3.cln.all %>% 
  mutate(., Ion.m = substr(Ion, 2, nchar(Ion)-1))

lipid.no.conc <- Exp3.all %>%
  filter(., !Class %in% IS_w_conc.Exp3$Class) #60 lipids have no IS for quantification

lipid.conc <- Exp3.all %>%
  filter(., Class %in% IS_w_conc.Exp3$Class) #354 lipids have IS for quantification

lipidmelt.df <- lipid.conc %>%   
  ungroup() %>% 
  select(., c(Class, LipidIon, Ion.m, QC_1:QC_3, XZ_83:XZ_104)) 

lipid.w.2xcal.chol <- bind_rows(lipidmelt.df, Xcal.chol.df)

lipidmelt <- lipid.w.2xcal.chol %>% ## 8900 rows
  pivot_longer(-c(LipidIon, Class, Ion.m), names_to = "Samples", values_to = "Intensity")


##append IS value of each sample with corresponding spike in concentration
##single point concentration calculation = [detected_lipid_intensity] /([internal_standard_intensity]/[concentration_of_internal_standard])
df.conc <- lipidmelt %>%
  mutate(., smpl = Samples) %>%
  mutate(., Class1 = Class) %>%
  mutate(., Ion1 = Ion.m) %>%
  group_by(Class, Ion.m, Samples) %>%
  group_modify(~ {
    .x %>%
      mutate(., IS_int = as.numeric(IS_w_conc.Exp3[IS_w_conc.Exp3$Class == unique(Class1) 
                                                    & IS_w_conc.Exp3$Ion == unique(Ion1), 
                                                    unique(smpl)])) %>%
      mutate(., conc_in_smpl = as.numeric(IS_w_conc.Exp3[IS_w_conc.Exp3$Class == unique(Class1) 
                                                          & IS_w_conc.Exp3$Ion == unique(Ion1), "Conc_in_Batch3"])) %>%
      mutate(., Endo_conc = Intensity/(IS_int/conc_in_smpl))
  })

### sanity check
s.n <- "XZ_88"
Cla.n <- "PE"
r.ind <- 38
a<- df.conc$Intensity[df.conc$Class == Cla.n & df.conc$Samples == s.n][r.ind]/df.conc$Endo_conc[df.conc$Class == Cla.n & df.conc$Samples == s.n][r.ind]
b <- as.numeric(IS_w_conc.Exp3[IS_w_conc.Exp3$Class == Cla.n, s.n])/IS_w_conc.Exp3$Conc_in_Batch3[IS_w_conc.Exp3$Class == Cla.n]
a == b

df.conc.wide <- df.conc %>%
  ungroup() %>% 
  select(., c("LipidIon", "Endo_conc", "Samples")) %>%
  pivot_wider(names_from = Samples, values_from = Endo_conc)

Exp3.lipid.w.conc <- df.conc.wide
save(Exp3.lipid.w.conc, file = paste0("./Output_data/Exp3_354_conc_lipidsw2endo.chol_356.Rdata")) 

Exp3.lipid.no.conc.c <- lipid.no.conc %>% 
  ungroup() %>% 
  select(., c(LipidIon, QC_1:QC_3, XZ_83:XZ_104))

Exp3_lpd.all <- bind_rows(Exp3.lipid.w.conc, Exp3.lipid.no.conc.c)#416 = 354+60+2
save(Exp3_lpd.all, file = paste0("./Output_data/Exp3_356_conc_60_noconc_lipids.Rdata"))


#### ========================================================#### ========================================================#### ========================================================
#### ========================================================#### ========================================================#### ========================================================
### O8-aNSC sample was accidentally diluted 4 times during sample preparation.
### To account for the difference, adjust concentration and intensity for O8-aNSC sample

lipid.w.conc.O8A.adj <- Exp3.lipid.w.conc %>% 
  mutate(., XZ_97 = 	XZ_97 *4)
save(lipid.w.conc.O8A.adj, file = paste0("./Output_data/Exp3_354_conc_lipidsw2endo.chol_356_O8A-adj.Rdata")) 

lipid.no.conc.c.O8A.adj <- Exp3.lipid.no.conc.c %>% 
  mutate(., XZ_97 = 	XZ_97 *4)

Exp3_lpd.all.O8A.adj <- bind_rows(lipid.w.conc.O8A.adj, lipid.no.conc.c.O8A.adj)#416 = 354+60+2
save(Exp3_lpd.all.O8A.adj, file = paste0("./Output_data/Exp3_356_conc_60_noconc_lipids_O8A-adj.Rdata"))
