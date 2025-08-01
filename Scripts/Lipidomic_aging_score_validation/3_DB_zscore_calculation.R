# Z score calculation on double bond composition features
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(stringi)

## ---------------------------------------------------------------------------------------------------------------------
# Load DB pct analysis on the dataset
load("./Output_Data/Ding.et.al.all.brain.DB.pct.Rdata")
load("./Output_Data/Meta_DB_signature.Rdata")

br.db.sgnt <- br.DB.PCT %>% 
  filter(Cla_DB %in% DB.hi.old)

#' 
#' #### Z score on DB PCT 
#' Z score by Cla_DB only, across all area of the brain
## ---------------------------------------------------------------------------------------------------------------------
z.DB.sgnt <- br.db.sgnt %>% 
  group_by(Cla_DB) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(DB_Pct))
  })

Z.sc.sgnt.DB <- z.DB.sgnt %>% 
  group_by(Samples) %>% 
  summarise(Mean.DB_Pct.Z = mean(zscore, na.rm = TRUE)) %>% 
  rowwise() %>% 
  mutate(., Age = str_split(Samples, ":")[[1]][5]) %>% 
  mutate(., Sex = str_split(Samples, ":")[[1]][6]) %>% 
  mutate(., Tissue = str_split(Samples, ":")[[1]][4])

Z.sc.sgnt.DB$Age <- factor(Z.sc.sgnt.DB$Age, levels = c("3 weeks", "16 weeks", "59 weeks", "92 weeks"))
save(Z.sc.sgnt.DB, file = "./Output_Data/Ding.et.all.DB.zscore.on.aging.signature.Rdata")

