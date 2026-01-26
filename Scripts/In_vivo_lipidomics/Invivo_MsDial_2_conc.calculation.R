### Concentration calculation

# In vivo sorted data
setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)
library(openxlsx)

load("./Output_Data/Invivo_IS_w_conc.Rdata")
Invivo.df <- read.csv("./Input_Data/Invivo_sorted_NSCs_integration_graded_120325.csv", stringsAsFactors = F)
Invivo.peak <- Invivo.df %>% 
  filter(Peak_shape == "Yes") #29

df.fix <- Invivo.peak %>% 
  mutate(Metabolite.name = ifelse(Alignment.ID == 3476, "ST 27:1;O", Metabolite.name)) %>% 
  mutate(Metabolite = ifelse(Alignment.ID == 3476, "ST 27:1;O", Metabolite)) 

Invivo.lpd.check <- df.fix %>% 
  rowwise() %>% 
  mutate(Class = str_split(Metabolite.name, " ")[[1]][1]) %>%
  mutate(Class = ifelse(Class == "ST", "Cholesterol", Class)) %>% 
  relocate(Class, .after = Metabolite.name) %>% 
  group_by(Metabolite.name) %>% 
  group_modify(~{
    .x %>% 
      mutate(Ion.check = ifelse(Adduct.type == Invivo.IS.Conc$Adduct.type[Invivo.IS.Conc$Class == Class], "T", "F"))  
  }) %>% 
  relocate(c(Class, Ion.check), .after = "Metabolite.name") %>% 
  filter(Ion.check == "T") ##all 29 lipids check out


## Append sample name to file name
label <- read.csv("./Input_Data/Lable_sample_ID.csv", stringsAsFactors = F)
label.df <- label %>% 
  filter(!Running_ID == "XY.18")
Invivo.endo.int <- Invivo.lpd.check %>% 
  ungroup() %>% 
  rename_with(~ str_remove(.x, "pos_")) %>% 
  select(c(Metabolite, matches("^XY."), Adduct.type)) %>% 
  rename_at(vars(matches(label.df$Running_ID)), ~label.df$Sorted.cell.samples) %>% 
  select(c(Metabolite, Adduct.type, matches("^Y|^O"))) %>% 
  pivot_longer(-c(Metabolite, Adduct.type), names_to = "Sample", values_to = "Intensity") %>% 
  rowwise() %>% 
  mutate(Class = str_split(Metabolite, " ")[[1]][1]) %>% 
  mutate(Class = ifelse(Class == "ST", "Cholesterol", Class))

##append IS value of each sample with corresponding spike in concentration
##single point concentration calculation = [detected_lipid_intensity] /([internal_standard_intensity]/[concentration_of_internal_standard])

df.conc <- Invivo.endo.int %>%
  mutate(., smpl = Sample) %>%
  mutate(., Class1 = Class) %>%
  mutate(., Ion1 = Adduct.type) %>%
  group_by(Class, Adduct.type, Sample) %>%
  group_modify(~ {
    .x %>%
      mutate(., IS_int = as.numeric(Invivo.IS.Conc[Invivo.IS.Conc$Class == unique(Class1) 
                                                   & Invivo.IS.Conc$Adduct.type == unique(Ion1), 
                                                   unique(smpl)])) %>%
      mutate(., conc_in_smpl = as.numeric(Invivo.IS.Conc[Invivo.IS.Conc$Class == unique(Class1) 
                                                         & Invivo.IS.Conc$Adduct.type  == unique(Ion1), "Conc_in_sample"])) %>%
      mutate(., Endo_conc = Intensity/(IS_int/conc_in_smpl))
  })


### sanity check
s.n <- "O4"
Cla.n <- "LPE"
r.ind <- 1
a<- df.conc$Intensity[df.conc$Class == Cla.n & df.conc$Sample == s.n][r.ind]/df.conc$Endo_conc[df.conc$Class == Cla.n & df.conc$Sample == s.n][r.ind]
b <- as.numeric(Invivo.IS.Conc[Invivo.IS.Conc$Class == Cla.n, s.n])/Invivo.IS.Conc$Conc_in_sample[Invivo.IS.Conc$Class == Cla.n]
round(a,5) == round(b,5)

df.conc.wide <- df.conc %>%
  ungroup() %>% 
  select(., c("Sample", "Metabolite", "Endo_conc")) %>%
  pivot_wider(names_from = Sample, values_from = Endo_conc)

Invivo.lipid.w.conc.MsD.29peak <- df.conc.wide
save(Invivo.lipid.w.conc.MsD.29peak, file = paste0("./Output_data/Invivo_29.lipids.conc.byMsD.Rdata"))
