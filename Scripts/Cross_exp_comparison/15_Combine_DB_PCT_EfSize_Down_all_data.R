##Decrease with Age
##Aging double bond composition effect size ranking between datasets

##Input: 2 invitro datasets + in vivo + GPMV, for Fig S5

##Steps:
#1. Import double bond composition from all data sets
#2. Rank all Class:DB that are decreased with age based on combined ranking from all 4 datasets
#3. Plot consistently down-regulated Class:DB in descending order based on effect size.
#4. Save Cla:DB list to create aging signature for later scripts
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")

##Primary Culture #1
load("./Output_Data/Ef_Size_DB_pct_InVitro.Rdata")
##Primary Culture #2
load("./Output_Data/Ef_Size_DB_PCT_Exp2_all_KO.Rdata")
#Subset to only include control sample
Exp2.ctrl <- Exp2.DB.es.g.allKO %>%
  filter(., KO == "N")
##In vivo
load("./Output_Data/Ef_Size_DB_pct_Invivo.Rdata")
##GPMV
load("./Output_Data/Ef_Size_DB_pct_GPMV.Rdata")

###Down with age####
LC.down <- down.ES.DB.func(LC.Invitro.DB.es.g, Cla_DB, es_g)
Invivo.down <- down.ES.DB.func(Invivo.DB.es.g, Cla_DB, es_g)
KO.ntg.down <- down.ES.DB.func(Exp2.ctrl, Cla_DB, es_g)
GPMV.down <- down.ES.DB.func(GPMV.DB.es.g, Cla_DB, es_g)

all.down <- list(LC.down, KO.ntg.down, Invivo.down, GPMV.down) %>% 
  reduce(full_join, by = "Cla_DB")
col.from <- c("Rank_pos.x", "Effect_Size.x", "Rank_Pct.x", "Rank_pos.y", "Effect_Size.y", "Rank_Pct.y",
              "Rank_pos.x.x", "Effect_Size.x.x", "Rank_Pct.x.x", "Rank_pos.y.y", "Effect_Size.y.y", "Rank_Pct.y.y")
col.to <- c("Rank_pos.In_vitro", "Effect_Size.In_vitro", "Rank_Pct.In_vitro",
            "Rank_pos.KO_In_vitro", "Effect_Size.KO_In_vitro", "Rank_Pct.KO_In_vitro",
            "Rank_pos.In_vivo", "Effect_Size.In_vivo", "Rank_Pct.In_vivo",
            "Rank_pos.GPMV", "Effect_Size.GPMV", "Rank_Pct.GPMV")

##calculate combined ranking in merged list, convert NA to 100, lower combined ranking means higher in the ranked list ####
all.down.cmb <- all.down
all.down.cmb[is.na(all.down.cmb)] = 100
all.down.cmb <- all.down.cmb %>%
  rename_at(vars(all_of(col.from)), function(x) col.to) %>%
  group_by(Cla_DB) %>%
  group_modify(~ {
    .x %>%
      mutate(., Cmb_Rank = sum(c(Rank_Pct.In_vitro, Rank_Pct.KO_In_vitro, Rank_Pct.In_vivo, Rank_Pct.GPMV), na.rm = T))
  }) %>%
  arrange(Cmb_Rank)

##get number of overlap from the merged list####
all.down.ovlp <- all.down %>%
  ungroup() %>%
  mutate(., no_overlap = rowSums(is.na(all.down))) %>%
  arrange(., no_overlap) %>%
  rename_at(vars(col.from), function(x) col.to)
all.down.ovlp$Cla_DB <- factor(all.down.ovlp$Cla_DB, levels = rev(all.down.cmb$Cla_DB))

overlap_order <- all.down.ovlp %>%
  select(., c("Cla_DB", "no_overlap"))

all.down.plot <- all.down.ovlp %>%
  select(., matches("Cla_DB|Rank_Pct")) %>%
  pivot_longer(-Cla_DB, names_to = "Data Sets", values_to = "Rank_Percentile") %>%
  mutate(., `Data Sets` = str_replace(`Data Sets`, "Rank_Pct.", ""))

all.down.plot <- left_join(all.down.plot, overlap_order, by = "Cla_DB")   %>% 
  mutate(., `Data Sets` = case_when(
  `Data Sets` == "In_vitro" ~ "Primary \nqNSC culture \n#1",
  `Data Sets` == "KO_In_vitro" ~ "Primary \nqNSC culture \n#2",
  `Data Sets` == "In_vivo" ~ "In vivo \nisolated \nqNSCs",
  `Data Sets` == "GPMV" ~ "GPMV",
))

all.down.plot$`Data Sets` <- factor(all.down.plot$`Data Sets`, levels = c("Primary \nqNSC culture \n#1", 
                                                                      "Primary \nqNSC culture \n#2", 
                                                                      "In vivo \nisolated \nqNSCs",
                                                                      "GPMV"))

##Subset to only keep consistent DB composition change in at least three out of four datasets####
all.down.plot.but1 <- all.down.plot %>% 
  filter(., no_overlap <=3) #68 entries in all datasets, 17 Class: DB combo
p <- ggplot(data = all.down.plot.but1,
            aes(x= `Data Sets`, y= Cla_DB, fill= Rank_Percentile)) 
p +  geom_tile(height = 1) +
  scale_fill_gradient2(low = "#26408B", high = "white", mid = "#A0D2E7", 
                       midpoint = 50)+
  labs(fill= "Rank")+
  labs(title = "Down DB in 3 or more datasets (w GPMV)", x = "", y = "Double Bond Abundance", color = "")+
  theme(text=element_text(size =13, face = "bold"))+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5),
        axis.text = element_text(size = 12, face = "plain", colour = "black"))
ggsave(paste0("./Figure_Panels/Fig_S5d_down.pdf"), width = 5, height = 5, useDingbats=FALSE)


### Check data sparsity and label annotate with "not detected"
DB.list <- factor(unique(all.down.plot.but1$Cla_DB))
DB.ls.dw <- unique(as.character(all.down.plot.but1$Cla_DB))
LC.check <- LC.Invitro.DB.es.g %>%
  filter(., Cla_DB %in% DB.ls.dw) %>%
  select(., c(Cla_DB, es_g)) %>%
  rename(., "exp1" = "es_g")

LC2.check <- Exp2.ctrl %>%
  filter(., Cla_DB %in% DB.ls.dw) %>%
  select(., c(Cla_DB, es_g)) %>%
  rename(., "exp2" = "es_g")

Invivo.check <- Invivo.DB.es.g %>%
  filter(., Cla_DB %in% DB.ls.dw) %>%
  select(., c(Cla_DB, es_g)) %>%
  rename(., "invivo" = "es_g")

GPMV.check <- GPMV.DB.es.g %>%
  filter(., Cla_DB %in% DB.ls.dw) %>%
  select(., c(Cla_DB, es_g)) %>%
  rename(., "gpmv" = "es_g")

check.all <- list(LC.check, LC2.check, Invivo.check,GPMV.check) %>% 
  reduce(full_join, by = "Cla_DB")

check.all$Cla_DB <- factor(check.all$Cla_DB, levels = rev(levels(DB.list)))
check.all <- check.all %>%
  arrange(.,Cla_DB)

write.csv(check.all, file = paste0("./Output_Data/Down_DB_Ef_Size_all_data_check_sparsity.csv"))
  
