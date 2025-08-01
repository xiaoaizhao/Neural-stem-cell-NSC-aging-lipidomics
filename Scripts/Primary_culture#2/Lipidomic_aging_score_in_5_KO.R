# Lipidomic aging score in all KOs from Primary culture #2

## --------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(rstatix)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")
## -------------------------------------------------------------------------------------------------------------------
load('Output_Data/Exp2_DB_PCT_all_samples.Rdata') # double bond composition in Primary Culture #2 Exp with KO
load("./Output_Data/Exp2_Norm_Impt_backtoraw_all693_lipids.Rdata") # individual lipids in Primary Culture #2 Exp with KO
load("./Output_Data/KO_LUT_paper.order2024-02-08.Rdata") # color scheme

## Import lipidomic aging signature ------
load("./Output_Data/Meta_Lipid_signature.Rdata")
load("./Output_Data/Meta_DB_signature.Rdata")

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
  rownames_to_column(., var = "LipidIon")

lpd.zs.df <- ether.rename(lpd.zs) %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") %>%
  group_by(LipidIon) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(Conc_Int))
  })

highold.mean.lpd.zs <- lpd.zs.df %>%
  filter(., LipidIon %in% conc.lpd.hi.old) 

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

Cmb.avg$KO.name <- factor(Cmb.avg$KO.name, levels = KO.LUT.paperorder$KO)
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
## result: all p-values are above 0.05, does not reject the hypothesis of normality. Will perform t-test statistical analysis

Cmb.avg$KO <- factor(Cmb.avg$KO, 
                      levels = c("N", 
                                 "E", 
                                 "M",
                                 "F", 
                                 "A", 
                                 "P"))

Cmb.avg$Age <- factor(Cmb.avg$Age, levels = c("Young", "Old"))



## Pair-wise t-test for old and young samples between KO and control conditions
## -------------------------------------------------------------------------------------------------------------------
stat.test <- Cmb.avg %>%
  group_by(Age) %>%
  t_test(Cmb.avg.lpd.db.mean.zs ~ KO, paired = TRUE, ref.group = "N") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "KO")
stat.test$p.adj <- round(stat.test$p.adj, digits = 4)


a<- ggplot(Cmb.avg, aes(x= KO, y= Cmb.avg.lpd.db.mean.zs, color = KO.name))
a+ geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.18, alpha=0.8, size=3, shape = 2)+
  scale_color_manual(values = KO.LUT.paperorder$Color)+
  scale_shape_manual(values =c(`Young` = 16, `Old` = 17)) +
  theme_classic()+
  theme(text=element_text(size = 13, face = "plain"),
        axis.text = element_text(colour = "black"))+
  labs(title = "Lipidomic aging score" , x = "", y = "Mean Z Score", color = "")+
  facet_wrap(vars(Age), nrow = 1) +
  stat_pvalue_manual(stat.test, label = "p.adj",
                     bracket.nudge.y = 0.1,
                     step.increase = 0.03, size = 3.5)+
  theme(legend.position= "none")

ggsave(filename = "./Figure_Panels/Fig.4e.pdf", width = 5, height = 5, useDingbats=FALSE)


## -------------------------------------------------------------------------------------------------------------------

Cmb.avg.med <- Cmb.avg %>% 
  group_by(Age_KO) %>% 
  summarise(., MedianZscore = median(Cmb.avg.lpd.db.mean.zs))

Cmb.avg.med

