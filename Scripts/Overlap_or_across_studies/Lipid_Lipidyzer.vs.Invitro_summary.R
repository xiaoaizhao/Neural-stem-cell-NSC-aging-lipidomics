
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(stringi)
library(tidyverse)

load("./Output_Data/Lipids.2outof3LC.invitro.Sig.features.Rdata")

load("./Output_Data/Ldz_Ef_Size_Age_qNSC.Rdata")# Lipidyzer
load("./Output_data/Lipidyzer_qNSC_for.LC-MS.ovlp.Rdata")
Ldz.ovlp <- Lipidyzer.Age.es.g %>% 
  filter(Lipid %in% Ldz.Qui.nodup$LipidIon) %>% 
  rowwise() %>% 
  mutate(ID_string = Ldz.Qui.nodup$ID_string[Ldz.Qui.nodup$LipidIon == Lipid])

CI95 <- function(df, mean, sem){
  df <- df %>% 
    rowwise() %>% 
    mutate(`Margin of Error` = 1.96 * {{sem}}) %>% 
    mutate(CI.lower = {{mean}} -`Margin of Error`) %>% 
    mutate(CI.upper = {{mean}} +`Margin of Error`)
}

c2.3.sig <- lpd.c2.in.3lc.df.sig %>% 
  select(LipidIon, MeanES) %>% 
  group_by(LipidIon) %>% 
  group_modify(~{
    .x %>% 
      summarise(SumES = unique(MeanES))
  }) %>% 
  rowwise() %>% 
  mutate(., ID_string = list(str_split(LipidIon, "\\(|_|\\)")[[1]][1:length(str_split(LipidIon, "\\(|_|\\)")[[1]])-1]))

Ldz.2of3 <- list(c2.3.sig, 
                  Ldz.ovlp
                  ) %>% 
  reduce(inner_join, by = c("ID_string")) %>%  #51 lipids overlap with In vivo
  rename("ES_summary" = "SumES") %>% 
  rename("ES_Lipidyzer" = "es_g")


Ldz.2of3.l <- Ldz.2of3 %>% 
  ungroup() %>% 
  select(ID_string, starts_with("ES")) %>% 
  pivot_longer(-ID_string, names_to = "Exp", values_to = "Effect size")

lpd.Ldz.2of3 <- Ldz.2of3.l%>% 
  group_by(ID_string) %>% 
  summarise(MeanES = mean(`Effect size`, na.rm = TRUE), SEM = sd(`Effect size`, na.rm = TRUE)/sqrt(sum(!is.na(`Effect size`)))) 

lpd.Ldz.2of3.CI.LoUp <- CI95(lpd.Ldz.2of3, MeanES, SEM)


c2of3all.p <- inner_join(Ldz.2of3.l, lpd.Ldz.2of3.CI.LoUp, by = "ID_string") %>% 
  mutate(Sig = ifelse( CI.lower > 0 & CI.upper > 0 | CI.lower < 0 & CI.upper < 0,
         "Significant", "Not significant")) %>% 
  group_by(ID_string) %>% 
  group_modify(~{
    .x %>% 
      mutate(star.pos = case_when(
        Sig == "Significant" & MeanES > 0 ~  max(`Effect size`) + 0.8, 
        Sig == "Significant" & MeanES < 0 ~ min(`Effect size`) - 0.8)
)
  }) %>% 
  mutate(Exp = case_when(
    grepl("summary", Exp) ~ "In vitro summary",
    grepl("Lipidyzer", Exp) ~ "Lipidyzer",
  )) %>% 
  rowwise() %>% 
  mutate(Lipid = paste(unlist(ID_string), collapse =  "_")) %>% 
  mutate(LipidIon = stri_replace_first(Lipid, "(", fixed = "_")) %>% 
  mutate(LipidIon = paste0(LipidIon, ")")) #66
  

c2of3all.p$Exp <- factor(c2of3all.p$Exp, levels = c( "In vitro summary", "Lipidyzer"))


a <-ggplot(c2of3all.p, aes(x = fct_reorder(LipidIon, MeanES), y = `Effect size`))
a+
  geom_point(aes(shape = Exp), colour = "grey39", alpha = 0.85, size = 3.5) +
  scale_shape_manual(values = c(13, 10)) + 
  geom_errorbar(aes(ymin = MeanES - SEM, ymax = MeanES + SEM), colour = "grey15", alpha = 0.75, width = 0.2) +
  stat_summary(aes(x=LipidIon,y=MeanES), fun=mean, geom = "point", size=3.5, shape=20, alpha = 0.75, colour = "grey15") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black", size = 8), axis.text.x = element_text(angle = 0, vjust = 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() + 
  labs(title = "Lipidyzer vs. 2 out of 3LC - significant lipids", x = "", y = "Effect size (Old vs. Young)") +
  geom_point(data = c2of3all.p %>% 
           filter(Sig == "Significant"), 
           aes(x = fct_reorder(LipidIon, MeanES), y = star.pos),
         pch=8, 
         size=1.2, stroke = 0.7, alpha = 0.75,
         colour="black") +
  theme(legend.position = "bottom")
ggsave(filename = "./Figure_Panels/EDFig.3h.pdf", width = 4, height = 6, useDingbats = FALSE)

