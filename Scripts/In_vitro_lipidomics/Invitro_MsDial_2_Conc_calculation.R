### Concentration calculation

setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)
library(stringi)

Invitro.new.peak <- read.csv("./Input_Data/Exp3_Validated_Reintegration_121125_color.code.csv", stringsAsFactors = F)

Invitro <- Invitro.new.peak %>% 
  filter(!Peak.color.code == "orange") #404 lipids
load("./Output_Data/Invitro.MS-Dial_IS.all.Rdata")

Invitro.lpd.check <- Invitro %>% 
  rowwise() %>% 
  mutate(Class = str_split(Metabolite.name, " ")[[1]][1]) %>% 
  mutate(Class = ifelse(Class == "ST", "Cholesterol", Class)) %>% 
  group_by(Metabolite.name) %>% 
  group_modify(~{
    .x %>% 
      mutate(Ion.check = ifelse(Adduct.type == IS.invitro.export.flt$Adduct.type[IS.invitro.export.flt$Class == Class], "T", "F"))  
  }) %>% 
  relocate(c(Class, Ion.check), .after = "Metabolite.name") %>% 
  filter(Ion.check == "T") ##all 404 lipids check out


## transform into long form
Lpd.invitro <- Invitro.lpd.check %>% 
  select(c(Metabolite.name, Metabolite, Adduct.type, matches("XZ|QC"))) %>% 
  pivot_longer(-c(Metabolite.name, Adduct.type, Metabolite), names_to = "Sample", values_to = "Intensity") %>% 
  rowwise() %>% 
  mutate(Class = str_split(Metabolite.name, " ")[[1]][1]) %>% 
  mutate(Class = ifelse(Class == "ST", "Cholesterol", Class)) %>% 
  filter(!grepl("PBS|QC", Sample))


##append IS value of each sample with corresponding spike in concentration
##single point concentration calculation = [detected_lipid_intensity] /([internal_standard_intensity]/[concentration_of_internal_standard])
load("./Output_data/IS_Invitro_w_conc_MsD.Rdata")
df.conc <- Lpd.invitro %>%
  mutate(., smpl = Sample) %>%
  mutate(., Class1 = Class) %>%
  mutate(., Ion1 = Adduct.type) %>%
  group_by(Class, Adduct.type, Sample) %>%
  group_modify(~ {
    .x %>%
      mutate(., IS_int = as.numeric(IS_w_conc.Invitro[IS_w_conc.Invitro$Class == unique(Class1) 
                                                   & IS_w_conc.Invitro$Adduct.type == unique(Ion1), 
                                                   unique(smpl)])) %>%
      mutate(., conc_in_smpl = as.numeric(IS_w_conc.Invitro[IS_w_conc.Invitro$Class == unique(Class1) 
                                                         & IS_w_conc.Invitro$Adduct.type  == unique(Ion1), "Conc_in_Batch3"])) %>%
      mutate(., Endo_conc = Intensity/(IS_int/conc_in_smpl))
  })

# 	Cholesterol IS was not detected in XZ_102, XZ_102 is Y2_aNSC-A
df.st <- df.conc %>% #Cholesterol IS was not detected in XZ_102, so the calculated endogenous IS is infinite. Will manually change this to 0
 mutate(Endo_conc = ifelse(Class == "Cholesterol" & smpl == "XZ_102", 0, Endo_conc))


### sanity check
s.n <- "XZ_100"
Cla.n <- "PC"
r.ind <- 8
a<- df.st$Intensity[df.st$Class == Cla.n & df.st$Sample == s.n][r.ind]/df.st$Endo_conc[df.st$Class == Cla.n & df.st$Sample == s.n][r.ind]
b <- as.numeric(IS_w_conc.Invitro[IS_w_conc.Invitro$Class == Cla.n, s.n])/IS_w_conc.Invitro$Conc_in_Batch3[IS_w_conc.Invitro$Class == Cla.n]
a == b

df.conc.wide <- df.st %>%
  ungroup() %>% 
  select(., c("Sample", "Metabolite.name", "Metabolite", "Endo_conc")) %>%
  pivot_wider(names_from = Sample, values_from = Endo_conc)

Invitro.lipid.w.conc.MsD <- df.conc.wide
save(Invitro.lipid.w.conc.MsD, file = paste0("./Output_data/Invitro_404.quant.lipids.byMsD.Rdata"))


##==== Attach sample name to each file ====
Smp.Key <- read.csv("./Input_Data/Batch3_sample_key_forR.csv", stringsAsFactors = F)
Smp.Key.e <- Smp.Key %>% 
  select(Sample.Name, Sample_ID) %>% 
  filter(!grepl("PBS", Sample_ID))

## Annotate sample ID with sample name
Invitro.lpd.conc <- Invitro.lipid.w.conc.MsD %>% 
  rename_at(vars(matches(Smp.Key.e$Sample_ID)), ~Smp.Key.e$Sample.Name)

## ==== Adjust concentration for O8-aNSC sample ==== 
### O8-aNSC sample was accidentally diluted 4 times during sample preparation.
### To account for the difference, adjust concentration and intensity for O8-aNSC sample

Invitro.lpd.conc.O8aNSC.adj <- Invitro.lpd.conc %>% 
  mutate(`O8_aNSC-A` = `O8_aNSC-A` *4  )

save(Invitro.lpd.conc.O8aNSC.adj, file = paste0("./Output_data/Invitro_404.quant.lipids.Smpl.Name.O8aNSC.adj.Rdata"))
