### Comparative lipid analysis between In vitro and In vitro Experiment #2

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)

source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Lpd.MsD.Age.ES_Exp2_all_KO.w.RT.Rdata") 
load("./Output_Data/Lpd.MsD.Invitro_Qui_Age_ES.w.RT.Rdata")
load("./Output_Data/Common.lipids_E2_E3.Rdata") #175 common lipid when unique identifier is stripped

E2 <- MsD.lpd.rmv.abc(Exp2.lpd.es.g.MsD.allKO.wRT) %>% 
  filter(KO == "N") %>% 
  select(LipidIon, Average.Rt.min., es_g) %>% 
  filter(LipidIon %in% common.E2.E3) #192

Invitro <- MsD.lpd.rmv.abc(Invitro.Q.lpd.ES.w.RT) %>% 
  select(LipidIon, Average.Rt.min., es_g) %>% 
  filter(LipidIon %in% common.E2.E3) #204

Invitro.E2.merge <- full_join(E2, Invitro, by = "LipidIon")  %>% #225
  rename_with(~ str_replace_all(., "\\.x", "Exp2")) %>% 
  rename_with(~ str_replace_all(., "\\.y", "Exp3")) 

Invitro.E2.tally <- Invitro.E2.merge %>% 
  group_by(LipidIon) %>% 
  tally() %>% 
  mutate(Unique = ifelse(n == 1, "Yes", "No")) %>% 
  group_split(Unique)

UniqueInvitroE2 <- Invitro.E2.merge %>% 
  filter(LipidIon %in% Invitro.E2.tally[[2]]$LipidIon) #134

## ==== function to test if all values are identical
are_all_identical <- function(x) {
  # Remove NA if you want to ignore them; keep if NA should break equality
  x <- x[!is.na(x)]
  length(unique(x)) <= 1
}


DupInvitroE2 <- Invitro.E2.merge %>% 
  filter(LipidIon %in% Invitro.E2.tally[[1]]$LipidIon) %>% #121
  mutate(DeltaRT = abs(Average.Rt.min.Exp3 - Average.Rt.min.Exp2)) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
    .x %>% 
      mutate(Keep_1 = ifelse(any(are_all_identical(Average.Rt.min.Exp2), are_all_identical(Average.Rt.min.Exp3)), "Yes", "No"))
  })

Dup.keep1 <- DupInvitroE2 %>% 
  filter(Keep_1 == "Yes") %>% #75
  group_by(LipidIon) %>% 
  group_modify( ~{
    .x %>% 
      filter(DeltaRT == min(DeltaRT))
  }) #37

## ==== For case of multi - multi duplicates matching, export and manualy sort in Excel ====
Dup.multi <- DupInvitroE2 %>% 
  filter(Keep_1 == "No") 

Dup.clean <- read.csv(file = "./Input_Data/Invitro.vs.Exp2.sorted.csv", stringsAsFactors = F)  
Dup.clean <- Dup.clean[,-1]

ES.invitro.all <- bind_rows(UniqueInvitroE2, Dup.keep1, Dup.clean) #175 total = 134 overlapping lipid + 37 with dups (cleaned based on minimum delta RT) + 4 lipids from manual sorting
save(ES.invitro.all, file = "./Output_Data/Overlapping.2.invitro.lpds.w.RT.n.ES.Rdata")

CI95 <- function(df, mean, sem){
  df <- df %>% 
    rowwise() %>% 
    mutate(`Margin of Error` = 1.96 * {{sem}}) %>% 
    mutate(CI.lower = {{mean}} -`Margin of Error`) %>% 
    mutate(CI.upper = {{mean}} +`Margin of Error`)
}

ES.sum <- ES.invitro.all %>% 
  select(LipidIon, es_gExp2, es_gExp3) %>% 
  pivot_longer(-LipidIon, names_to = "Exp", values_to = "Effect size") %>% 
  group_by(LipidIon) %>% 
  summarise(MeanES = mean(`Effect size`, na.rm = TRUE), SEM = sd(`Effect size`, na.rm = TRUE)/sqrt(sum(!is.na(`Effect size`)))) 

c2.in.2lc.lpd.CI.LoUp <- CI95(ES.sum, MeanES, SEM)

E2.p <-ES.invitro.all %>% 
  select(LipidIon, es_gExp2) %>% 
  rename("es_g" = "es_gExp2") %>% 
  mutate(Exp = "In vitro Experiment #2")

Invitro.p <- ES.invitro.all %>% 
  select(LipidIon, es_gExp3) %>% 
  rename("es_g" = "es_gExp3") %>% 
  mutate(Exp = "In vitro")


c2.in.2lc.df <- bind_rows(E2.p, Invitro.p)
p.c2.in.2lc.df <- left_join(c2.in.2lc.df, c2.in.2lc.lpd.CI.LoUp, by = "LipidIon") %>% 
   mutate(Sig = ifelse( CI.lower > 0 & CI.upper > 0 | CI.lower < 0 & CI.upper < 0,
         "Significant", "Not significant"))

sig.2of2lc <- p.c2.in.2lc.df %>% 
   filter(Sig == "Significant") #160

non.sig <- p.c2.in.2lc.df %>% 
  filter(Sig == "Not significant") #190
 
sig.2of2lc$Exp <- factor(sig.2of2lc$Exp, levels = c( "In vitro", "In vitro Experiment #2"))


## ==== Calculate p value based on CI ====
df.2of2.pval.CI.t <- pval.from.CI(sig.2of2lc, 0.05, 0) ## 160 entries

lpd.c2.in.2lc.df.sig.t.final <- df.2of2.pval.CI.t %>% 
  filter(padj < 0.05) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
    .x %>% 
      mutate(star.pos = case_when(
        MeanES > 0 ~  max(es_g) + 0.8, 
        MeanES < 0 ~ min(es_g) - 0.8)
      ) 
  }) #160 entires - this means all p value significant lipids are also FDR significant

save(lpd.c2.in.2lc.df.sig.t.final, file = "./Output_Data/Lipids.2of2LC.invitro.FDRSig.features.final.Rdata") 

a <-ggplot(lpd.c2.in.2lc.df.sig.t.final, aes(x = fct_reorder(LipidIon, MeanES), y = es_g))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 2, fill = "gold") +
  scale_shape_manual(values = c(17,2, 10)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM), colour = "grey15", alpha = 0.75, width = 0.2) +
  stat_summary(aes(x=LipidIon,y=MeanES), fun=mean, geom = "point", size=2, shape=20, alpha = 0.75, colour = "grey15") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 7.5), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() + 
  labs(title = "Significant in 2 in vitro studies", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = lpd.c2.in.2lc.df.sig.t.final %>% 
               filter(Sig == "Significant"), 
             aes(x = fct_reorder(LipidIon, MeanES), y = star.pos),
             pch=8, 
             size=0.8, stroke = 0.7, alpha = 0.75,
             colour="black") +
  theme(legend.position = "bottom")
ggsave(filename = paste0("./Figure_Panels/Fig.1e.pdf"), width =4, height = 8.5, useDingbats = FALSE)


