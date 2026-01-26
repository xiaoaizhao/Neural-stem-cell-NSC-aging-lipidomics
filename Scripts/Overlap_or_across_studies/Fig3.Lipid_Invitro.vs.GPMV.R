
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)

source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Lpd.MsD.GPMV_Age_ES.Rdata")
load("./Output_Data/Lpd.MsD.Invitro_Qui_Age_ES.w.RT.Rdata")
load("./Output_Data/Common.lipids_GPMV_Invitro.Rdata")

common.no.TG.DG <- common.GPMV.Invitro[!grepl("^TG|^DG", common.GPMV.Invitro)] 

GPMV <- MsD.lpd.rmv.abc(GPMV.lpd.ES.w.RT) %>%
  select(LipidIon, Average.Rt.min., es_g) %>%
  filter(LipidIon %in% common.no.TG.DG) #183

Invitro <- MsD.lpd.rmv.abc(Invitro.Q.lpd.ES.w.RT) %>%
  select(LipidIon, Average.Rt.min., es_g) %>%
  filter(LipidIon %in% common.no.TG.DG) #198

GPMV.Invitro.merge <- full_join(GPMV, Invitro, by = "LipidIon")  %>% #242
  rename_with(~ str_replace_all(., "\\.x", "GPMV")) %>% 
  rename_with(~ str_replace_all(., "\\.y", "Exp3")) 

GPMVInvitro.tally <- GPMV.Invitro.merge %>% 
  group_by(LipidIon) %>% 
  tally() %>% 
  mutate(Unique = ifelse(n == 1, "Yes", "No")) %>% 
  group_split(Unique)

UniqueGPMVInvitro <- GPMV.Invitro.merge %>% 
  filter(LipidIon %in% GPMVInvitro.tally[[2]]$LipidIon) #128

## ==== function to test if all values are identical
are_all_identical <- function(x) {
  # Remove NA if you want to ignore them; keep if NA should break equality
  x <- x[!is.na(x)]
  length(unique(x)) <= 1
}


DupGPMVInvitro <- GPMV.Invitro.merge %>% 
  filter(LipidIon %in% GPMVInvitro.tally[[1]]$LipidIon) %>% #94
  mutate(DeltaRT = abs(Average.Rt.min.Exp3 - Average.Rt.min.GPMV)) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
    .x %>% 
      mutate(Keep_1 = ifelse(any(are_all_identical(Average.Rt.min.GPMV), are_all_identical(Average.Rt.min.Exp3)), "Yes", "No"))
  })

Dup.keep1 <- DupGPMVInvitro %>% 
  filter(Keep_1 == "Yes") %>% #62
  group_by(LipidIon) %>% 
  group_modify( ~{
    .x %>% 
      filter(DeltaRT == min(DeltaRT))
  }) #31

## ==== For case of multi - multi duplicates matching, export and manualy sort in Excel ====
Dup.multi <- DupGPMVInvitro %>% 
  filter(Keep_1 == "No") #32


Dup.clean <- read.csv(file = "./Input_Data/GPMV.vs.Invitro_sorted.csv", stringsAsFactors = F)  
Dup.clean <- Dup.clean[,-1]

ES.GPMVInvitro.all <- bind_rows(UniqueGPMVInvitro, Dup.keep1, Dup.clean) #167 total = 128 overlapping lipid + 31 with dups (cleaned based on minimum delta RT) + 8 lipids from manual sorting
save(ES.GPMVInvitro.all, file = "./Output_Data/Overlapping.GPMVInvitro.lpds.w.RT.n.ES.Rdata")

CI95 <- function(df, mean, sem){
  df <- df %>% 
    rowwise() %>% 
    mutate(`Margin of Error` = 1.96 * {{sem}}) %>% 
    mutate(CI.lower = {{mean}} -`Margin of Error`) %>% 
    mutate(CI.upper = {{mean}} +`Margin of Error`)
}

ES.sum <- ES.GPMVInvitro.all %>% 
  select(LipidIon, es_gGPMV, es_gExp3) %>% 
  pivot_longer(-LipidIon, names_to = "Exp", values_to = "Effect size") %>% 
  group_by(LipidIon) %>% 
  summarise(MeanES = mean(`Effect size`, na.rm = TRUE), SEM = sd(`Effect size`, na.rm = TRUE)/sqrt(sum(!is.na(`Effect size`)))) 

c2.in.2lc.lpd.CI.LoUp <- CI95(ES.sum, MeanES, SEM)

GPMV.p <-ES.GPMVInvitro.all %>% 
  select(LipidIon, es_gGPMV) %>% 
  rename("es_g" = "es_gGPMV") %>% 
  mutate(Exp = "GPMV")

Invitro.p <- ES.GPMVInvitro.all %>% 
  select(LipidIon, es_gExp3) %>% 
  rename("es_g" = "es_gExp3") %>% 
  mutate(Exp = "In vitro")


c2.in.2lc.df <- bind_rows(GPMV.p, Invitro.p)
p.c2.in.2lc.df <- left_join(c2.in.2lc.df, c2.in.2lc.lpd.CI.LoUp, by = "LipidIon") %>% 
  mutate(Sig = ifelse( CI.lower > 0 & CI.upper > 0 | CI.lower < 0 & CI.upper < 0,
                       "Significant", "Not significant"))

sig.2of2lc <- p.c2.in.2lc.df %>% 
  filter(Sig == "Significant") #130

non.sig <- p.c2.in.2lc.df %>% 
  filter(Sig == "Not significant") #204

sig.2of2lc$Exp <- factor(sig.2of2lc$Exp, levels = c( "GPMV", "In vitro"))


## ==== Calculate p value based on CI ====
df.2of2.pval.CI.t <- pval.from.CI(sig.2of2lc, 0.05, 0) ## 130 entries

lpd.c2.in.2lc.df.sig.t.final <- df.2of2.pval.CI.t %>% 
  filter(padj < 0.05) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
    .x %>% 
      mutate(star.pos = case_when(
        MeanES > 0 ~  max(es_g) + 0.8, 
        MeanES < 0 ~ min(es_g) - 0.8)
      ) 
  }) #130 entires - this means all p value significant lipids are also FDR significant

GPMV.Invitro.lpd.fdr.sig <- lpd.c2.in.2lc.df.sig.t.final
save(GPMV.Invitro.lpd.fdr.sig, file = "./Output_Data/Lipids.GPMV.Invitro.FDRSig.features.Rdata") #65 features


### ==== plot all significant changes for extended figure ==== ####
a <-ggplot(GPMV.Invitro.lpd.fdr.sig, aes(x = fct_reorder(LipidIon, MeanES), y = es_g))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 2.3) +
  scale_shape_manual(values = c(15, 17)) +
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM), colour = "grey15", alpha = 0.75, width = 0.2) +
  stat_summary(aes(x=LipidIon,y=MeanES), fun=mean, geom = "point", size=2.6, shape=20, alpha = 0.75, colour = "grey15") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 5.5), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(title = "GPMV vs. In vitro", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = GPMV.Invitro.lpd.fdr.sig %>%
               filter(Sig == "Significant"),
             aes(x = fct_reorder(LipidIon, MeanES), y = star.pos),
             pch=8,
             size=1.2, stroke = 0.7, alpha = 0.75,
             colour="black") +
  theme(legend.position = "bottom")
ggsave(filename = paste0("./Figure_Panels/EDFig.8d.pdf"), width = 4, height = 9, useDingbats = FALSE)

### ==== Plot the top 30% biggest effect size in main

top30.GPMV <- GPMV.Invitro.lpd.fdr.sig %>% 
  group_by(LipidIon) %>% 
  summarise(MeanES = unique(MeanES)) %>% 
  mutate(percentile = percent_rank(abs(MeanES))) %>% 
  filter(percentile > 0.7) #20 lipids

top30.GPMV.p <- GPMV.Invitro.lpd.fdr.sig %>% 
  filter(LipidIon %in% top30.GPMV$LipidIon) #40 entries

a <-ggplot(top30.GPMV.p, aes(x = fct_reorder(LipidIon, MeanES), y = es_g))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 3) +
  scale_shape_manual(values = c(15, 17)) +
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM), colour = "grey15", alpha = 0.75, width = 0.2) +
  stat_summary(aes(x=LipidIon,y=MeanES), fun=mean, geom = "point", size=3, shape=20, alpha = 0.75, colour = "grey15") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 5.5), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(title = "Top 30% GPMV vs. In vitro", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = top30.GPMV.p %>%
               filter(Sig == "Significant"),
             aes(x = fct_reorder(LipidIon, MeanES), y = star.pos),
             pch=8,
             size=1.2, stroke = 0.7, alpha = 0.75,
             colour="black") +
  theme(legend.position = "bottom")
ggsave(filename = paste0("./Figure_Panels/Fig.3b.pdf"), width = 5, height = 5, useDingbats = FALSE)
