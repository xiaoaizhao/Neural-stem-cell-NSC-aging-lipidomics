### -------------------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(RColorBrewer)
library(ggrepel)
library(stringi)


## ----------------------------------------------------------------------------------------------------------------------------------
noslash <- function(x, na.rm = FALSE){str_replace_all(x, "/", "_")}


## ----------------------------------------------------------------------------------------------------------------------------------
df.all <- read.csv("./Input_Data/Brain_cell_lipidomics_Fitzner.csv", stringsAsFactors = F)
b.cla <- unique(df.all$class)
b.cla


## ----------------------------------------------------------------------------------------------------------------------------------
load("./Output_Data/Invitro.MsD.404.log2.norm.impt.no.O8aNSC.Rdata")
source("./Scripts/Function_scripts/Pre-processing_functions.R")


InvitroAQ <- 2^Invitro.no.O8A.MsD %>% 
  rownames_to_column(var = "LipidIon")

my.dt <- MsD.lpd.rmv.abc(InvitroAQ) %>% 
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>% 
  mutate(Class = ifelse(LipidIon == "Cholesterol", "Cholesterol", Class)) %>% 
  mutate(., SideChain = substr(LipidIon,
                               str_locate(LipidIon, "\\(")+1, 
                               str_locate(LipidIon, "\\)")-1))

my.cla <- unique(my.dt$Class)
my.cla


## ----------------------------------------------------------------------------------------------------------------------------------
c.cla <- intersect(my.cla, b.cla)
c.cla


## ----------------------------------------------------------------------------------------------------------------------------------
LC.TG <- InvitroAQ %>%
  filter(grepl("^TG", LipidIon)) %>%  #67
  mutate(Class = "TG") %>% 
  rowwise() %>% 
  mutate(SideChain = substr(LipidIon, 4, str_locate(LipidIon, "\\|")-1)) %>% 
  mutate(., ID_string = list(c(Class, SideChain)))



## ----------------------------------------------------------------------------------------------------------------------------------
LC.Cer.SM <- InvitroAQ %>%
  filter(grepl("^Cer|^SM", LipidIon))  %>%#27
  mutate(Class = substr(LipidIon, 1, str_locate(LipidIon, " ")-1)) %>% 
  rowwise() %>% 
  mutate(SideChain = substr(LipidIon, str_locate(LipidIon, " ")+1, str_locate(LipidIon, "\\;")-1)) %>% 
  mutate(., ID_string = list(c(Class, SideChain)))



## ----------------------------------------------------------------------------------------------------------------------------------

LC.LPC.LPE <- my.dt %>% #59
  filter(grepl("^LPE|^LPC", LipidIon)) %>% 
  mutate(., ID_string = list(c(Class, SideChain)))
  


## ----------------------------------------------------------------------------------------------------------------------------------
ls <- c("DG", "PC", "PE", "PG", "PI", "PS")
LC.2chain <- my.dt %>% 
  ungroup() %>%
  filter(Class %in% ls) %>%  #250
  rowwise() %>% 
  mutate(SepChain = list(str_split(SideChain, "_")[[1]]))

LC.2chain.df <- LC.2chain %>% 
  ungroup() %>% 
  mutate(ID_string = map2(Class, SepChain, ~ c(.x, .y))) %>% 
  relocate(ID_string, .after = "LipidIon") %>% 
  select(-SepChain)


## ----------------------------------------------------------------------------------------------------------------------------------
LC.chol <- my.dt %>% 
  ungroup() %>% 
  filter(LipidIon == "Cholesterol") %>% 
  rowwise() %>% 
  mutate(ID_string = map2(Class, "", ~ c(.x, .y)))


## ----------------------------------------------------------------------------------------------------------------------------------
LC.all <- bind_rows(LC.Cer.SM,
                    LC.2chain.df,
                    LC.TG,
                    LC.LPC.LPE,
                    LC.chol) %>% 
  arrange(LipidIon)#403 lipids (404 lipids total. - cholesterol)


## ----------------------------------------------------------------------------------------------------------------------------------
LC.AQ <- LC.all %>% 
  ungroup() %>% 
  relocate(ID_string, .after = "LipidIon") %>%
  select(., -LipidIon)

save(LC.AQ, file = "./Output_Data/LC_Invitro_for_admixture.Rdata")


## ----------------------------------------------------------------------------------------------------------------------------------
BR.smpl <- colnames(df.all)[str_detect(colnames(df.all), "Neurons|Oligo")]


## ----------------------------------------------------------------------------------------------------------------------------------
BR.ether <- df.all %>% 
  filter(., grepl("-", class))

BR.ether.1c <- BR.ether %>% 
  filter(., !grepl("/", feature)) %>% 
  mutate(., LipidIon = feature) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
  mutate(., Class = substr(feature, 1, str_locate(feature, " ")-1)) %>% 
  mutate(., SideChain = substr(feature, 
                               str_locate(feature, "O"),
                               str_locate(feature, ";")-1)) %>% 
  mutate(., ID_string = list(c(Class, SideChain)))
  }) %>% 
  ungroup() %>% 
  relocate(ID_string, .after = "feature") %>% 
  select(., c("ID_string","feature", BR.smpl))

BR.ether.2c <- BR.ether %>% 
  filter(., grepl("/", feature)) %>% 
  mutate(., LipidIon = feature) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
  mutate(., Class = substr(feature, 1, str_locate(feature, " ")-1)) %>% 
  mutate(., SideChain1 = substr(feature, 
                               str_locate(feature, "O"),
                               str_locate(feature, ";")-1)) %>% 
  mutate(., SideChain2 = substr(feature, 
                           str_locate(feature, "/")+1,
                           stri_locate_last(feature, fixed = ";")-1)) %>% 
  mutate(., ID_string = list(c(Class, SideChain1, SideChain2)))
  }) %>% 
  ungroup() %>% 
  relocate(ID_string, .after = "feature") %>% 
  select(., c("ID_string","feature", BR.smpl))
  


## ----------------------------------------------------------------------------------------------------------------------------------
BR.1c <- df.all %>% 
  filter(., !grepl("_|/", feature)) %>% 
  filter(., !feature %in% BR.ether.1c$feature) %>% 
  mutate(., LipidIon = feature) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
  mutate(., Class = substr(feature, 1, str_locate(feature, " ")-1)) %>% 
  mutate(., SideChain = substr(feature, 
                               str_locate(feature, " ")+1,
                               str_locate(feature, ";")-1)) #%>% 
  # mutate(., ID_string = list(c(Class, SideChain)))
}) %>% 
  mutate(., Class = ifelse(grepl("Chol", feature), "Cholesterol", Class)) %>% 
  mutate(., Class = ifelse(grepl("TAG", feature), "TG", Class)) %>% 
  mutate(., SideChain = ifelse(grepl("Chol", feature), "", SideChain)) %>% 
  mutate(., ID_string = list(c(Class, SideChain))) %>% 
  ungroup() %>% 
  relocate(ID_string, .after = "feature") %>% 
  select(., c("ID_string","feature", BR.smpl))



## ----------------------------------------------------------------------------------------------------------------------------------
BR.2c <- df.all %>% 
  filter(., grepl("_", feature)) %>% 
  filter(., !feature %in% BR.ether.2c$feature) %>% 
  mutate(., LipidIon = feature) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
  mutate(., Class = substr(feature, 1, str_locate(feature, " ")-1)) %>% 
  mutate(., SideChain1 = substr(feature, 
                               str_locate(feature, " ")+1,
                               str_locate(feature, ";")-1)) %>% 
  mutate(., SideChain2 = substr(feature, 
                           str_locate(feature, "_")+1,
                           stri_locate_last(feature, fixed = ";")-1)) #%>%
  # mutate(., ID_string = list(c(Class, SideChain1, SideChain2)))
}) %>% 
  mutate(., Class = ifelse(grepl("DAG", feature), "DG", Class)) %>% 
  mutate(., ID_string = list(c(Class, SideChain1, SideChain2))) %>% 
  ungroup() %>% 
  relocate(ID_string, .after = "feature") %>% 
  select(., c("ID_string","feature", BR.smpl))


## ----------------------------------------------------------------------------------------------------------------------------------
BR.all <- bind_rows(BR.1c, BR.2c, BR.ether.1c, BR.ether.2c) %>% 
  arrange(feature)

df.all <- df.all %>% 
  arrange(feature)
all.equal(BR.all$feature, df.all$feature)


## ----------------------------------------------------------------------------------------------------------------------------------
BR.Neu.Oligo <- BR.all %>% 
  ungroup() %>% 
  # filter(rowSums(is.na(BR.all)) < 9) %>% 
  select(., -feature)

save(BR.Neu.Oligo, file = "./Output_Data/Neuron_Oligo_for_admixture.Rdata")


## ----------------------------------------------------------------------------------------------------------------------------------

mix.cell <- full_join(LC.AQ, BR.Neu.Oligo, by = "ID_string") #868 404 (lc in vitro) + 682 (neuron/oligo) - 868 = 218 overlapping lipids

