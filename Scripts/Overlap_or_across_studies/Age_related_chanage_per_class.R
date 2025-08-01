
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(stringi)
source("./Scripts/Function_scripts/Effect_size_functions.R")

## ==== Function for LC-MS/MS datasets====
lpd.age.change <- function(df, dataset.name, fig.num.) {
df1 <- df %>% 
  pivot_longer(-LipidIon, names_to = "Samples", values_to = "Concentration") %>% 
  mutate(Age = ifelse(grepl("^Y", Samples), "Young", "Old")) %>% 
  group_by(LipidIon, Age) %>% 
  summarise(mean.conc = mean(Concentration))

df2 <- df1 %>% 
  pivot_wider(names_from = Age, values_from = mean.conc) %>% 
  mutate(Dlt.age = Old - Young)

df3 <- df2 %>% 
  ungroup() %>% 
  summarise(age.sum = sum(abs(Dlt.age)))

df4 <- df2 %>%
  mutate(age.sum.noDirction = df3$age.sum) %>% 
  mutate(age.pct = abs(Dlt.age/age.sum.noDirction) *100) %>% 
  mutate(Class = ifelse(grepl("Cholesterol", LipidIon),
                        "Cholesterol",
                        substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)))

df5 <- df4 %>% 
  group_by(Class) %>% 
  summarise(Class.sum.per.age = sum(age.pct))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")

mycolors <- lut.march25 %>% 
  filter(lipid_cat %in%df5$Class)

df5$Class <- factor(df5$Class, levels = c(mycolors$lipid_cat))

df5 <- df5 %>% 
    arrange(Class)

b <- ggplot(df5, aes(x="", y=Class.sum.per.age))+
  geom_bar(aes(fill=Class),
  position="stack", stat="identity")+
    labs(title  =  dataset.name)+
    ylab("mol% change with age")+
    xlab("")+
    scale_fill_manual(values = as.character(mycolors$Clr_list.1.18.)) +
    theme(legend.position= "right")+
    guides(fill=guide_legend(title="Class")) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"))

b

ggsave(filename = paste0("./Figure_Panels/",fig.num., ".pdf"), width = 3, height = 5, useDingbats = FALSE)
}  


## ==== Function for lipidyzer dataset====

ldz.lpd.age.change <- function(df, dataset.name, fig.num.) {
df1 <- df %>% 
  rownames_to_column(var = "LipidIon") %>% 
  pivot_longer(-LipidIon, names_to = "Samples", values_to = "Concentration") %>% 
  mutate(Age = ifelse(grepl("^Y", Samples), "Young", "Old")) %>% 
  group_by(LipidIon, Age) %>% 
  summarise(mean.conc = mean(Concentration))

df2 <- df1 %>% 
  pivot_wider(names_from = Age, values_from = mean.conc) %>% 
  mutate(Dlt.age = Old - Young)

df3 <- df2 %>% 
  ungroup() %>% 
  summarise(age.sum = sum(abs(Dlt.age)))

df4 <- df2 %>%
  mutate(age.sum.noDirction = df3$age.sum) %>% 
  mutate(age.pct = abs(Dlt.age/age.sum.noDirction) *100) %>% 
  mutate(Class = substr(LipidIon, 1, stri_locate_first(LipidIon, regex = "\\.")-1)) %>% 
  mutate(Class = ifelse(grepl("TAG", Class), "TG", Class)) %>% 
  mutate(Class = ifelse(grepl("DAG", Class), "DG", Class)) %>% 
  mutate(Class = ifelse(grepl("CER", Class), "Cer", Class)) %>% 
  mutate(Class = ifelse(Class == "CE", "ChE", Class)) %>% 
  filter(!Class == "FFA")

df5 <- df4 %>% 
  group_by(Class) %>% 
  summarise(Class.sum.per.age = sum(age.pct)) %>% 
  arrange(desc(Class.sum.per.age))

load("./Output_Data/Class_col_list_paper.order_031725.Rdata")

mycolors <- lut.march25 %>% 
  filter(lipid_cat %in%df5$Class)

df5$Class <- factor(df5$Class, levels = c(mycolors$lipid_cat))

df5 <- df5 %>% 
    arrange(Class)

b <- ggplot(df5, aes(x="", y=Class.sum.per.age))+
  geom_bar(aes(fill=Class),
  position="stack", stat="identity")+
  labs(title  =  dataset.name)+
  ylab("mol% change with age")+
  xlab("")+
  scale_fill_manual(values = as.character(mycolors$Clr_list.1.18.)) +
  theme(legend.position= "right")+
  guides(fill=guide_legend(title="Class")) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"))

b

ggsave(filename = paste0("./Figure_Panels/",fig.num., ".pdf"), width = 3, height = 5, useDingbats = FALSE)
}  

## ==== Primary culture #2====
load("./Output_Data/Exp2_Norm_Impt_backtoraw_conc612_lipids.Rdata")
KO.ctrl <- raw_conc.exp2 %>% 
  select(contains("_N")) %>% 
  rownames_to_column(var = "LipidIon")

KO.ctrl <- ether.rename(KO.ctrl)

lpd.age.change(KO.ctrl, "Primary_culture_#2", "EDFig.1g")

## ==== Primary culture #3====
load("./Output_Data/M2PM.contrl.samples.conc.lipids.Rdata")
lpd.age.change(M2PM.lpd.ctrl, "Primary_culture_#3", "Fig.1d")

## ==== In vivo====
load("./Output_Data/Invivo_Norm_Impt_RAW.conc_121_lipids.Rdata")
lpc.to.rmv <- c("LPC(14:0)", "LPC(16:0e)", "LPC(16:1e)", "LPC(16:2e)", "LPC(17:1)", "LPC(18:1e)", "LPC(18:3e)", "LPC(20:1)", "LPC(20:5)")

Invivo.df <- Invivo.Impt.norm.conc.raw %>% 
  mutate(LipidIon = substr(LipidIon, 1, str_locate(LipidIon, "\\)"))) %>% 
  filter(!LipidIon %in% lpc.to.rmv)

Invivo<- ether.rename(Invivo.df) 

lpd.age.change(Invivo, "In_vivo", "EDFig.5e")

## ==== GPMV====
load(file = "./Output_Data/GPMV_Norm_Impt_log2_conc_483_lipids.Rdata")
GPMV.df <- 2^GPMV_Impt_norm_conc_all %>% 
  rownames_to_column(., var = "LipidIon")

GPMV.df <- ether.rename(GPMV.df)

lpd.age.change(GPMV.df, "GPMV", "EDFig.8d")

## ==== Lipidyzer ====
load("./Output_data/Ldz_511_lipids_backtoraw.conc_MedConc_norm_imputed.Rdata")

ldz.lpd.age.change(ldz.raw_conc, "Lipidyzer", "EDFig.3g")
# sig.pct = length(q.es.Ldz$LipidIon)/ length(Lipidyzer.Age.es.g$Lipid) *100
# sig.pct #[1] 15.06849

