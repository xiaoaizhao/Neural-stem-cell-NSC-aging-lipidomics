
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(patchwork)
source(file= "./Scripts/Function_scripts/Pre-processing_functions.R")
source(file= "./Scripts/Function_scripts/Effect_size_functions.R")

load("./Output_Data/Lpd.MsD.Invitro_Qui_Age_ES.w.RT.Rdata")
load("./Output_Data/Lpd.MsD.Age.ES_Exp2_all_KO.w.RT.Rdata")
load("./Output_Data/Lpd.MsD.GPMV_Age_ES.Rdata")
load("./Output_Data/Lpd.MsD.Invivo_Age_ES.w.RT.Rdata")

# summary effect size from meta analysis
load("./Output_Data/Lipid.summary.from.meta.analysis.Rdata") 
# features that are significantly higher in old
load("./Output_Data/Meta_Lipid_signature.Rdata")


##==== Forest plot on lipid features ====
E2.meta <- MsD.lpd.rmv.abc(Exp2.lpd.es.g.MsD.allKO.wRT) %>% 
  filter(LipidIon %in% conc.lpd.hi.old & KO == "N")  %>% 
  select(LipidIon, es_g, se_g) %>% 
  mutate(Exp = "In vitro Experiment #2")
  
Invitro.meta <- MsD.lpd.rmv.abc(Invitro.Q.lpd.ES.w.RT) %>% 
  select(LipidIon, es_g, se_g) %>% 
  mutate(Exp = "In vitro") %>% 
  filter(LipidIon %in% conc.lpd.hi.old) 

GPMV.meta <- MsD.lpd.rmv.abc(GPMV.lpd.ES.w.RT) %>% 
  select(LipidIon, es_g, se_g) %>% 
  mutate(Exp = "GPMV") %>% 
  filter(LipidIon %in% conc.lpd.hi.old)

Invivo.meta <- MsD.lpd.rmv.abc(Invivo.lpd.ES.w.RT) %>% 
  select(LipidIon, es_g, se_g) %>% 
  mutate(Exp = "In vivo") %>% 
  filter(LipidIon %in% conc.lpd.hi.old)


meta.l <- Lipid.summary %>% 
  rename(., 'es_g' = 'summary') %>%
  rename(., 'se_g' = 'se_summary') %>%
  dplyr::select(LipidIon, es_g, se_g) %>% 
  mutate(Exp = "Meta_analysis") %>% 
  filter(LipidIon %in% conc.lpd.hi.old)

Lpd.meta.all <- bind_rows(E2.meta, Invitro.meta,
                          Invivo.meta, GPMV.meta, meta.l)

Lpd.meta.all.plot <- left_join(Lpd.meta.all, meta.l, by = "LipidIon") %>% 
  rename("Exp_es" = "es_g.x") %>% 
  rename("Meta_es" = "es_g.y") %>% 
  mutate(Cat = ifelse(Exp.x == "Meta_analysis", "Summary_ES", "Experiments"))

Lpd.meta.all.plot$Exp.x <- factor(Lpd.meta.all.plot$Exp.x, levels = c( "In vitro", "In vitro Experiment #2", "In vivo", "GPMV", "Meta_analysis"))

df.n <- Lpd.meta.all.plot %>% 
  group_by(Exp.x, LipidIon) %>% 
  tally() %>% 
  mutate(Exp_lpd = paste0(Exp.x, "_", LipidIon)) %>% 
  filter(n>1)

dup.df <- Lpd.meta.all.plot %>% 
  mutate(Exp_lpd = paste0(Exp.x, "_", LipidIon)) %>% 
  filter(Exp_lpd %in% df.n$Exp_lpd) %>% 
  group_by(Exp_lpd) %>% 
  group_modify(~{
    .x %>% 
      filter(Exp_es == max(Exp_es))
  })

uni.df <- Lpd.meta.all.plot %>% 
  mutate(Exp_lpd = paste0(Exp.x, "_", LipidIon)) %>% 
  filter(!Exp_lpd %in% df.n$Exp_lpd) 

lpd.plot.df <- bind_rows(uni.df, dup.df)

lpd.ftr <-ggplot(lpd.plot.df, aes(x = fct_reorder(LipidIon, Meta_es), y = Exp_es))+
  geom_errorbar(aes(ymin = Meta_es - se_g.y, ymax = Meta_es + se_g.y), width = 0.13, colour = "magenta3", alpha = 0.85) +
  geom_point(aes( shape = Exp.x, color = Cat, size = Cat), alpha = 0.75, stroke = 0.4) +
  scale_shape_manual(values = c(17,2, 25, 15, 20)) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c( "grey39", "magenta3")) +
  scale_size_manual(values = c( 2.5, 4.5)) +
  coord_flip() + 
  labs(title = "Lipid aging signature", x = "", y = "Effect size in qNSCs (Old vs. Young)") +
  theme(legend.position = "right")
lpd.ftr

lipid.n <- lpd.plot.df %>% 
  group_by(LipidIon, Exp.x) %>% 
  tally() %>% 
  arrange(desc(n))

##==== Forest plot on DB features ====
# rm(list = ls())
load('Output_Data/DBPct.MsD.Invitro_Qui_Age_ES.Rdata')
load('Output_Data/DB.MsD.Age.ES_Exp2_all_KO.Rdata')
load('Output_Data/DBPct.MsD.Invivo_Age_ES.Rdata')
load('Output_Data/DBPct.MsD.GPMV_Age_ES.Rdata')

load("./Output_Data/DB.summary.from.meta.analysis.Rdata")
load("./Output_Data/Meta_DB_signature.Rdata")


E2.meta <- Exp2.MsD.DB.es.g.allKO %>% 
  filter(KO == "N") %>% 
  select(Cla_DB, es_g, se_g) %>% 
  mutate(Exp = "In vitro Experiment #2") %>% 
  filter(Cla_DB %in% conc.DB.hi.old)

Invitro.meta <- Invitro.Qui.MsD.CONC.DB.es.g %>% 
  select(Cla_DB, es_g, se_g) %>% 
  mutate(Exp = "In vitro") %>% 
  filter(Cla_DB %in% conc.DB.hi.old)

GPMV.meta <- GPMV.DB.AG.ES %>% 
  select(Cla_DB, es_g, se_g) %>% 
  mutate(Exp = "GPMV") %>% 
  filter(Cla_DB %in% conc.DB.hi.old)

Invivo.meta <- Invivo.DB.AG.ES %>% 
  select(Cla_DB, es_g, se_g) %>% 
  mutate(Exp = "In vivo") %>% 
  filter(Cla_DB %in% conc.DB.hi.old)

meta.l <- DB.summary %>% 
  rename(., 'es_g' = 'summary') %>%
  rename(., 'se_g' = 'se_summary') %>%
  dplyr::select(Cla_DB, es_g, se_g) %>% 
  mutate(Exp = "Meta_analysis") %>% 
  filter(Cla_DB %in% conc.DB.hi.old)

DB.meta.all <- bind_rows(E2.meta, Invitro.meta,
                          Invivo.meta, GPMV.meta, meta.l) %>%
  filter(!is.nan(es_g)) 

DB.meta.all.plot <- left_join(DB.meta.all, meta.l, by = "Cla_DB") %>% 
  rename("Exp_es" = "es_g.x") %>% 
  rename("Meta_es" = "es_g.y") %>% 
  mutate(Cat = ifelse(Exp.x == "Meta_analysis", "Summary_ES", "Experiments"))

DB.meta.all.plot$Exp.x <- factor(DB.meta.all.plot$Exp.x, levels = c("In vitro", "In vitro Experiment #2", "In vivo", "GPMV", "Meta_analysis"))


db.ftr <-ggplot(DB.meta.all.plot, aes(x = fct_reorder(Cla_DB, Meta_es), y = Exp_es)) +
  geom_errorbar(aes(ymin = Meta_es - se_g.y, ymax = Meta_es + se_g.y), width = 0.08, colour = "magenta3", alpha = 0.85) +
  geom_point(aes( shape = Exp.x, color = Cat, size = Cat), alpha = 0.75, stroke = 0.4) +
  scale_shape_manual(values = c(17, 2, 15, 20)) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c( "grey39", "magenta3")) +
  scale_size_manual(values = c( 3, 4.5)) +
  coord_flip() + 
  labs(title = "DB aging signature", x = "", y = "Effect size in qNSCs (Old vs. Young)") +
  theme(legend.position = "bottom")
db.ftr

lpd.ftr + db.ftr 
ggsave(filename = paste0("./Figure_Panels/Fig.4d.pdf"), width = 10, height = 8, useDingbats = FALSE)



