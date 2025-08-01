
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(patchwork)
source(file= "./Scripts/Function_scripts/Effect_size_functions.R")

load("./Output_Data/Exp3_Qui_Age_ES.conc.lipids.Rdata")
load("./Output_Data/Ef_Size_CONC.Lipid_Exp2_all_KO.Rdata")
load("./Output_Data/Ef_Size_CONC.Lipid_InVivo.Rdata")
load("./Output_Data/Ef_Size_Conc.Lipid_GPMV.Rdata")

# summary effect size from meta analysis
load("./Output_Data/Lipid.summary.from.meta.analysis.Rdata") 
# features that are significantly higher in old
load("./Output_Data/Meta_Lipid_signature.Rdata")

##==== Forest plot on lipid features ====
KO.meta <- ether.rename(Exp2.CONC.lpd.es.g.allKO) %>% 
  filter(LipidIon %in% conc.lpd.hi.old & KO == "N")  %>% 
  select(LipidIon, es_g, se_g) %>% 
  mutate(Exp = "Primary culture #2")
  
E3.meta <- E3.Q.conc.AG.ES %>% 
  select(LipidIon, es_g, se_g) %>% 
  mutate(Exp = "Primary culture #3") %>% 
  filter(LipidIon %in% conc.lpd.hi.old) 

GPMV.meta <- ether.rename(GPMV.conc.lpd.es.g) %>% 
  select(LipidIon, es_g, se_g) %>% 
  mutate(Exp = "GPMV") %>% 
  filter(LipidIon %in% conc.lpd.hi.old)

Invivo.meta <- ether.rename(Invivo.CONC.lpd.es.g) %>% 
  select(LipidIon, es_g, se_g) %>% 
  mutate(Exp = "In vivo") %>% 
  filter(LipidIon %in% conc.lpd.hi.old)


meta.l <- Lipid.summary %>% 
  rename(., 'es_g' = 'summary') %>%
  rename(., 'se_g' = 'se_summary') %>%
  dplyr::select(LipidIon, es_g, se_g) %>% 
  mutate(Exp = "Meta_analysis") %>% 
  filter(LipidIon %in% conc.lpd.hi.old)

Lpd.meta.all <- bind_rows(KO.meta, E3.meta,
                          Invivo.meta, GPMV.meta, meta.l)

Lpd.meta.all.plot <- left_join(Lpd.meta.all, meta.l, by = "LipidIon") %>% 
  rename("Exp_es" = "es_g.x") %>% 
  rename("Meta_es" = "es_g.y") %>% 
  mutate(Cat = ifelse(Exp.x == "Meta_analysis", "Summary_ES", "Experiments"))

Lpd.meta.all.plot$Exp.x <- factor(Lpd.meta.all.plot$Exp.x, levels = c( "Primary culture #2", "Primary culture #3", "In vivo", "GPMV", "Meta_analysis"))

lpd.ftr <-ggplot(Lpd.meta.all.plot, aes(x = fct_reorder(LipidIon, Meta_es), y = Exp_es))+
  geom_errorbar(aes(ymin = Meta_es - se_g.y, ymax = Meta_es + se_g.y), width = 0.2, colour = "magenta3", alpha = 0.85) +
  geom_point(aes( shape = Exp.x, color = Cat, size = Cat), alpha = 0.75, fill = "grey80", stroke = 0.8) +
  scale_shape_manual(values = c(2, 17, 25, 15, 20)) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c( "grey39", "magenta3")) +
  scale_size_manual(values = c( 1.5, 3.5)) +
  coord_flip() + 
  labs(title = "4 sep data Conc Lipid aging features", x = "", y = "Effect size in qNSCs (Old vs. Young)") +
  theme(legend.position = "none")
lpd.ftr


##==== Forest plot on DB features ====
load('./Output_Data/Ef_Size_CONC.DB_pct_Exp3_Qui_aging.Rdata')
load('./Output_Data/Ef_Size_DB_PCT_Exp2_all_KO.Rdata')
load('./Output_Data/Ef_Size_DB_pct_Invivo.Rdata')
load('./Output_Data/Ef_Size_DB_pct_GPMV.Rdata')

load("./Output_Data/DB.summary.from.meta.analysis.Rdata")
load("./Output_Data/Meta_DB_signature.Rdata") # 5 double bond features


E2.meta <- Exp2.DB.es.g.allKO %>% 
  filter(KO == "N") %>% 
  select(Cla_DB, es_g, se_g) %>% 
  mutate(Exp = "Primary culture #2") %>% 
  filter(Cla_DB %in% DB.hi.old)

E3.meta <- Exp3.Qui.CONC.DB.es.g %>% 
  select(Cla_DB, es_g, se_g) %>% 
  mutate(Exp = "Primary culture #3") %>% 
  filter(Cla_DB %in% DB.hi.old)

GPMV.meta <- GPMV.DB.es.g %>% 
  select(Cla_DB, es_g, se_g) %>% 
  mutate(Exp = "GPMV") %>% 
  filter(Cla_DB %in% DB.hi.old)

Invivo.meta <- Invivo.DB.es.g %>% 
  select(Cla_DB, es_g, se_g) %>% 
  mutate(Exp = "In vivo") %>% 
  filter(Cla_DB %in% DB.hi.old)

meta.l <- DB.summary %>% 
  rename(., 'es_g' = 'summary') %>%
  rename(., 'se_g' = 'se_summary') %>%
  dplyr::select(Cla_DB, es_g, se_g) %>% 
  mutate(Exp = "Meta_analysis") %>% 
  filter(Cla_DB %in% DB.hi.old)

DB.meta.all <- bind_rows(E2.meta, E3.meta,
                          Invivo.meta, GPMV.meta, meta.l) %>%
  filter(!is.nan(es_g)) 

DB.meta.all.plot <- left_join(DB.meta.all, meta.l, by = "Cla_DB") %>% 
  rename("Exp_es" = "es_g.x") %>% 
  rename("Meta_es" = "es_g.y") %>% 
  mutate(Cat = ifelse(Exp.x == "Meta_analysis", "Summary_ES", "Experiments"))

DB.meta.all.plot$Exp.x <- factor(DB.meta.all.plot$Exp.x, levels = c("Primary culture #2", "Primary culture #3", "In vivo", "GPMV", "Meta_analysis"))

DB.meta.all.plot.noSPHP <- DB.meta.all.plot %>% 
  filter(!Cla_DB == "SPHP:1")

db.ftr <-ggplot(DB.meta.all.plot.noSPHP, aes(x = fct_reorder(Cla_DB, Meta_es), y = Exp_es)) + 
  geom_errorbar(aes(ymin = Meta_es - se_g.y, ymax = Meta_es + se_g.y), width = 0.04, colour = "magenta3", alpha = 0.75) +
  geom_point(aes( shape = Exp.x, color = Cat, size = Cat), alpha = 0.75, fill = "grey80", stroke = 0.8) +
  scale_shape_manual(values = c(2, 17, 25, 15, 20)) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c( "grey39", "magenta3")) +
  scale_size_manual(values = c( 1.5, 3.5)) +
  coord_flip() + 
  labs(title = "4 sep data conc/noconc DB aging features", x = "", y = "Effect size in qNSCs (Old vs. Young)")
db.ftr

lpd.ftr + db.ftr 
ggsave(filename = "./Figure_Panels/Fig.4d.pdf", width = 10, height = 8, useDingbats = FALSE)
