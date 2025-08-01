## Double bond composition analysis in order to compute the lipidomic aging score based on side chain composition, in addition to lipid species
## ------------------------------------------------------------------
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(stringi)

#' 
#' ### Features that goes up with age
## ---------------------------------------------------------------------------------------------------------------------
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Ding.et.al.lipid.data.for.analysis.Rdata")

#' 
## ---------------------------------------------------------------------------------------------------------------------
DB <- br.lpd.fmt %>% 
  rowwise() %>% 
  mutate(SideChain = list(substr(unlist(mt.name), 
                                 str_locate(unlist(mt.name), "\\(")+1, 
                                 str_locate(unlist(mt.name), "\\)")-1))) %>% 
  mutate(Class = unique(substr(unlist(mt.name), 1, str_locate(unlist(mt.name), "\\(")-1)))
  

#' 
#' Tally side chain
## ---------------------------------------------------------------------------------------------------------------------
all.br.DB <- db.tally.list(DB, Conc, Samples)
all.br.DB_by_class <- dplyr::bind_rows(all.br.DB)

#' 
#' Class Sum and DB_Pct
## ---------------------------------------------------------------------------------------------------------------------
br.cla.sum <- DB %>% 
  group_by(Samples, Class) %>% 
  summarise(ClassSum = sum(Conc))

br.DB.PCT<- left_join(all.br.DB_by_class, br.cla.sum, by = c("Samples", "Class")) %>%
  mutate(., DB_Pct = Sum_DB/ClassSum) %>% 
  mutate(Cla_DB = paste0(Class, DB_num))

save(br.DB.PCT, file = "./Output_Data/Ding.et.al.all.brain.DB.pct.Rdata")

