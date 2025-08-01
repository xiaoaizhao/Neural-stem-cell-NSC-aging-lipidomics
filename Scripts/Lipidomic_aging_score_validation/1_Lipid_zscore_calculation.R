# Z score calculation on lipid features

setwd(rstudioapi::getActiveProject())

rm(list = ls())
library(tidyverse)
library(stringi)

load("./Output_Data/Ding.et.al.lipid.data.for.analysis.Rdata")

load("./Output_Data/Meta_Lipid_signature.Rdata")

#' 
#' Z score lipids - by individual species only across all areas of the brain
## ---------------------------------------------------------------------------------------------------------------------
Z.lpd.br <- br.lpd.fmt %>% 
  group_by(mt.name) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(Conc))
  })

Z.sgnt.lpd <- Z.lpd.br %>% 
  rowwise() %>% 
  mutate(Ft.lpd.check = any(unlist(mt.name %in% conc.lpd.hi.old)))

Z.sc.sgnt.lpd <- Z.sgnt.lpd %>% 
  filter(Ft.lpd.check == "TRUE") %>% 
  group_by(Samples) %>% 
  summarise(Mean.Lpd.Z = mean(zscore, na.rm = TRUE)) %>% 
  rowwise() %>% 
  mutate(., Age = str_split(Samples, ":")[[1]][5]) %>% 
  mutate(., Sex = str_split(Samples, ":")[[1]][6]) %>% 
  mutate(., Tissue = str_split(Samples, ":")[[1]][4])

save(Z.sc.sgnt.lpd, file = "./Output_Data/Ding.et.all.Lipid.zscore.on.aging.signature.Rdata")

