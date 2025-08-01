#' Organize data from Fiehn group aging brain lipidomics
##------------------------------------------------------------------
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(stringi)


#' 
#' Matadata of all samples
## ---------------------------------------------------------------------------------------------------------------------
# mt.df <- as.data.frame(t(read.csv("./Input_Data/Metadata_Fiehn.csv", stringsAsFactors = F)))
mt.df <- read.csv("./Input_Data/Metadata_Fiehn.csv", stringsAsFactors = F, header = F)
concat.func <- function(x){
  paste(x, collapse = ":")
}
mt.n <- apply(mt.df, 2, concat.func)
mt.un <- unname(mt.n)[-1]

#' 
#' Raw data
## ---------------------------------------------------------------------------------------------------------------------
br.all <- read.csv("./Input_Data/Fiehn_brain_aging_raw_data.csv", stringsAsFactors = F)

br.all.nm <- br.all %>% 
  rename_at(vars(starts_with("Concentration..unit..ng.g")), ~ mt.un) %>% 
  select(-c("X.", "Quantification.method"))

#' 
#' Organize data
## ---------------------------------------------------------------------------------------------------------------------
br.conc <- br.all.nm %>% 
  select(-Annotation) %>% 
  mutate_if(is.character, as.numeric)
br.df <- bind_cols(br.all.nm$Annotation, br.conc) %>% 
  rename("Metabolite" = "...1")  %>% 
  pivot_longer(-Metabolite, names_to = "Samples", values_to = "Conc") %>% 
  rowwise() %>% 
  mutate(., Age = str_split(Samples, ":")[[1]][5]) %>% 
  mutate(., Sex = str_split(Samples, ":")[[1]][6]) %>% 
  mutate(., Tissue = str_split(Samples, ":")[[1]][4])
  
unique(br.df$Tissue) #10 areas of the brain total
# "Basal ganglia"   "Cerebellum"      "Cerebral cortex" "Hippocampus"     "Hypothalamus"    "Medulla"         "Midbrain"       "Olfactory bulb"  "Pons"            "Thalamus"  

#' 
#' Re-format metabolite name, add class and side chain columns
#' ### Features that goes up with age
#' without subsetting tissue, same the matrix for DB analysis
## ---------------------------------------------------------------------------------------------------------------------
br.lpd <- br.df %>% 
  mutate(anno = str_split(Metabolite, "or ")) %>% 
  rowwise() %>% 
  mutate(anno1= list(str_replace_all(anno, "; ", ")"))) %>% 
  mutate(anno1= list(str_replace_all(anno1, ";", ")"))) %>% 
  mutate(anno1= list(str_replace_all(anno1, " ", "("))) 

br.phospho.lpd <- br.lpd %>% 
  rowwise() %>% 
  filter(any(grepl("^PC|^PS|^PE|^PI|^PG", anno1))) %>% 
  filter(any(unlist(str_detect(anno1, "_")))) %>%  
  mutate(mt.name = list(substr(unlist(anno1), stri_locate_last(unlist(anno1), regex = "PC|PS|PE|PI|PG"), nchar(unlist(anno1)))))

br.LPE <- br.lpd %>% 
  filter(any(grepl("^LPE", anno1))) %>% 
  rowwise() %>% 
  mutate(., anno = list(paste0(anno, ")"))) %>% 
  mutate(anno1= list(str_replace_all(anno, "; ", ""))) %>% 
  mutate(anno1= list(str_replace_all(anno1, " ", "("))) %>% 
  mutate(mt.name = list(anno1))

br.lpd.fmt <- bind_rows(br.phospho.lpd, br.LPE) 
save(br.lpd.fmt, file = "./Output_Data/Ding.et.al.lipid.data.for.analysis.Rdata")

