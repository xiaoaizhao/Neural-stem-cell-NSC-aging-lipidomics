## Calculate lipidomic aging score and compare across samples in Primary Culture #2 (i.e. KO experiment)
## Aging lipidomic score is calculated as the mean z score across aging features 
## Will do this in 
# 1) all features combined (DB composition + individual lipid)

## --------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(RColorBrewer)
library(rstatix)

setwd(rstudioapi::getActiveProject())
## -------------------------------------------------------------------------------------------------------------------
load('Output_Data/Exp2_DB_PCT_all_samples.Rdata') # effect size matrix on double bond composition in Primary Culture #2 Exp with KO
load("./Output_Data/Exp2_Norm_Impt_backtoraw_all693_lipids.Rdata") # effect size matrix on individual lipid in Primary Culture #2 Exp with KO
load("./Output_Data/KO_LUT.Rdata") # Color LUT for each KO condition
load("./Output_Data/Meta_DB_signature_hi_in_Old_4_studies.Rdata") # Aging double bond composition features from meta-analysis
load("./Output_Data/Meta_Lipid_signature_hi_in_Old_4_studies.Rdata") # Aging individual lipid features from meta-analysis

## Get z score for each aging double bond composition features 
## -------------------------------------------------------------------------------------------------------------------
DB.zs <- Exp2_DB %>%
  mutate(., Cla_DB = paste0(Class, DB_num)) %>%
  group_by(Cla_DB) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(DB_Pct))
  })

highold.mean.DB.zs <- DB.zs %>%
  filter(., Cla_DB %in% DB.hi.old) %>%
  group_by(Sample) %>%
  summarise(., mean_zs = mean(zscore)) %>%
  mutate(., Culture = substr(Sample, 1, str_locate(Sample, "_")-1)) %>%
  mutate(., KO = substr(Sample, str_locate(Sample, "_")+1, nchar(Sample))) %>%
  mutate(., Age = ifelse(grepl("Y", Culture), "Young", "Old")) %>% 
  rename(., "DB_mean_zs" = "mean_zs") 

## -------------------------------------------------------------------------------------------------------------------
lpd.zs <- raw_int.exp2 %>%
  rownames_to_column(., var = "LipidIon") %>%
  mutate(., LipidIon = str_replace_all(LipidIon, "/", "_")) %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") %>%
  group_by(LipidIon) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(Conc_Int))
  })

highold.mean.lpd.zs <- lpd.zs %>%
  filter(., LipidIon %in% lpd.hi.old) 

highold.mean.lpd.zs <- highold.mean.lpd.zs %>% 
  group_by(Sample) %>%
  summarise(., mean_zs = mean(zscore)) %>%
  mutate(., Culture = substr(Sample, 1, str_locate(Sample, "_")-1)) %>%
  mutate(., KO = substr(Sample, str_locate(Sample, "_")+1, nchar(Sample))) %>%
  mutate(., Age = ifelse(grepl("Y", Culture), "Young", "Old")) %>% 
  rename(., "Lpd_mean_zs" = "mean_zs") 

## -------------------------------------------------------------------------------------------------------------------
Cmb <- left_join(highold.mean.DB.zs, highold.mean.lpd.zs, by = c("Culture", "Sample", "KO", "Age")) 
  
Cmb.avg <- Cmb %>% 
  group_by(., Sample) %>% 
  group_modify(~{
    .x %>% 
      mutate(., Cmb.avg.lpd.db.mean.zs = mean(c(DB_mean_zs, Lpd_mean_zs)))
  }) %>% 
  mutate(., KO.name = case_when(
    KO == "N" ~ "Control",
    KO == "E" ~ "Elovl5",
    KO == "M" ~ "Mboat2",
    KO == "A" ~ "Agpat3",
    KO == "F" ~ "Fads2",
    KO == "P" ~ "Pla2g4e",
  )) %>%
  mutate(., Age_KO = paste0(Age, "_", KO))

Cmb.avg$Age <- factor(Cmb.avg$Age, levels = c("Young", "Old"))

Cmb.avg$KO.name <- factor(Cmb.avg$KO.name, levels = KO.LUT$KO)


## Normality test to determine if it's appropriate to use parametric statistical testing
## -------------------------------------------------------------------------------------------------------------------
## Test between control and each KO condition

norm.test <- list()
for (x in levels(Cmb.avg$KO.name)[2:6]){
  
d <- with(Cmb.avg, 
          Cmb.avg.lpd.db.mean.zs[KO.name == "Control"] - Cmb.avg.lpd.db.mean.zs[KO.name == x])

norm.test[[x]] <- shapiro.test(d)$p.value
}

norm.df <- bind_rows(norm.test)
print(norm.df)
## result: all p-values are above 0.05, does not reject the hypothesis of normality. Will perform t-test to test KO effect in old samples


## Subset to only include old samples and young control samples
## -------------------------------------------------------------------------------------------------------------------
Cmb.avg <- Cmb.avg %>% 
  filter(., Age == "Old" |
           Age == "Young" & KO == "N")
Cmb.avg$Age_KO <- factor(Cmb.avg$Age_KO, levels = c("Young_N", "Old_N", "Old_E", "Old_M","Old_F", "Old_A", "Old_P"))


old.cmb <- Cmb.avg %>% 
  group_by(Sample) %>% 
  filter(., Age == "Old") 

YO_N <- Cmb.avg %>% 
  group_by(Sample) %>% 
  filter(., KO == "N") 

## One-way ANOVA shows significant difference between KO conditions
## -------------------------------------------------------------------------------------------------------------------
paired.aov <-aov(Cmb.avg.lpd.db.mean.zs ~ Age_KO +  Error(Culture), data = old.cmb)
summary(paired.aov)

## Pair-wise t-test for old samples between KO condition
## Un-paired t-test between young and old control samples
## -------------------------------------------------------------------------------------------------------------------
stat.test <- old.cmb %>%
  group_by(Age) %>%
  t_test(Cmb.avg.lpd.db.mean.zs ~ Age_KO, paired = TRUE, ref.group = "Old_N") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "Age_KO")

stat.YO<- YO_N %>%
  group_by(KO) %>%
  t_test(Cmb.avg.lpd.db.mean.zs ~ Age_KO, paired = FALSE, ref.group = "Young_N") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.YO <- stat.YO %>% add_xy_position(x = "Age_KO")

stat.all <- bind_rows(stat.test, stat.YO)

a<- ggplot(Cmb.avg, aes(x= Age_KO, y= Cmb.avg.lpd.db.mean.zs, color = KO.name, shape = Age))
a+ geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.18, alpha=0.8, size=3)+
  scale_color_manual(values = KO.LUT$Color)+
  scale_shape_manual(values =c(`Young` = 16, `Old` = 17)) +
  theme_classic()+
  theme(text=element_text(size = 13, face = "plain"),
        axis.text = element_text(colour = "black"))+
  labs(title = "Lipid + DB composition aging signature combined" , x = "", y = "Mean Z Score", color = "")+
  stat_pvalue_manual(stat.all, label = "p.adj",
                     bracket.nudge.y = 0.1,
                     step.increase = 0.03, size = 3.5)

ggsave(filename = "./Figure_Panels/Fig_4d.pdf", width = 5, height = 5, useDingbats=FALSE)


## -------------------------------------------------------------------------------------------------------------------

Cmb.avg.med <- Cmb.avg %>% 
  group_by(Age_KO) %>% 
  summarise(., MedianZscore = median(Cmb.avg.lpd.db.mean.zs))

Cmb.avg.med
