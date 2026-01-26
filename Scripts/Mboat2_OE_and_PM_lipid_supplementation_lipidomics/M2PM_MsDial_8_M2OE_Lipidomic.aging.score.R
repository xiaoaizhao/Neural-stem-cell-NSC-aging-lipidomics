### Lipidomic aging score in Mboat2 overexpression and plasma membrane lipid supplementation
## --------------------------------------------------------------------
rm(list=ls())
library(rstatix)
library(tidyverse)
library(ggpubr)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")
source("./Scripts/Function_scripts/Pre-processing_functions.R")
## -------------------------------------------------------------------------------------------------------------------
load("./Output_Data/KO_LUT_paper.order2024-02-08.Rdata") # color scheme

## Import lipidomic aging signature ------
load("./Output_Data/Meta_Lipid_signature.Rdata")
load("./Output_Data/Meta_DB_signature.Rdata")

## Mboat2 OE samples 
load("./Output_Data/Mboat2_OE_DB_PCT.Rdata")
load("./Output_Data/Mboat2_OE_LIPID.Rdata")


## Get z score for each aging double bond composition features 
## -------------------------------------------------------------------------------------------------------------------
DB.zs <- M2OE.DB %>%
  mutate(., Cla_DB = paste0(Class, DB_num)) %>%
  group_by(Cla_DB) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(DB_Pct))
  }) %>% #976
  filter(!is.nan(zscore[,1])) #944

highold.mean.DB.zs <- DB.zs %>%
  filter(., Cla_DB %in% conc.DB.hi.old) %>% 
  group_by(Sample) %>%
  summarise(., mean_zs = mean(zscore, na.rm = TRUE)) %>%
  mutate(Condition = case_when(
    grepl("_Mb2_OE", Sample) ~ "Mboat2 OE",
    grepl("_EGFP", Sample) ~ "Control"
  )) %>% 
  mutate(., Age = ifelse(grepl("^Y", Sample), "Young", "Old")) %>% 
  rename(., "DB_mean_zs" = "mean_zs") %>% 
  rename("Sample.Name" = "Sample")

## -------------------------------------------------------------------------------------------------------------------
lpd.zs.n <-M2OE.Lipid %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample.Name", values_to = "Conc_Int") %>%
  group_by(LipidIon) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(Conc_Int))
  }) %>% 
  ungroup()

lpd.zs <- MsD.lpd.rmv.abc(lpd.zs.n)

highold.mean.lpd.zs <- lpd.zs %>%
  filter(., LipidIon %in% conc.lpd.hi.old) 

highold.mean.lpd.zs <- highold.mean.lpd.zs %>% 
  group_by(Sample.Name) %>%
  summarise(., mean_zs = mean(zscore, na.rm = TRUE)) %>%
  mutate(Condition = case_when(
    grepl("_Mb2_OE", Sample.Name) ~ "Mboat2 OE",
    grepl("_EGFP", Sample.Name) ~ "Control"
  )) %>% 
  mutate(., Age = ifelse(grepl("^Y", Sample.Name), "Young", "Old")) %>% 
  rename(., "Lpd_mean_zs" = "mean_zs") 

## -------------------------------------------------------------------------------------------------------------------
Cmb <- left_join(highold.mean.DB.zs, highold.mean.lpd.zs, by = c("Sample.Name", "Condition", "Age")) 

Cmb.avg <- Cmb %>% 
  group_by(., Sample.Name) %>% 
  group_modify(~{
    .x %>% 
      mutate(., Cmb.avg.lpd.db.mean.zs = mean(c(DB_mean_zs, Lpd_mean_zs*0.6)))
  }) 

Cmb.avg$Condition <- factor(Cmb.avg$Condition, levels = c("Control", "Mboat2 OE"))
Cmb.avg$Age <- factor(Cmb.avg$Age, levels = c("Young", "Old"))
## Normality test to determine if it's appropriate to use parametric statistical testing
## -------------------------------------------------------------------------------------------------------------------
## Test between control and treatment in each age group

norm.test <- list()
for (x in levels(Cmb.avg$Age)){
  
  d <- with(Cmb.avg, 
            Cmb.avg.lpd.db.mean.zs[Condition == "Control"] - Cmb.avg.lpd.db.mean.zs[Condition == "Mboat2 OE"])
  
  norm.test[[x]] <- shapiro.test(d)$p.value
}

norm.df <- bind_rows(norm.test)
print(norm.df)
## result: all p-values are above 0.05, does not reject the hypothesis of normality. Will perform t-test for 


## Pair-wise t-test for old and young samples between OE and control conditions
## -------------------------------------------------------------------------------------------------------------------
stat.test <- Cmb.avg %>%
  group_by(Age) %>%
  t_test(Cmb.avg.lpd.db.mean.zs ~ Condition, paired = TRUE, ref.group = "Control") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "Condition")

M2.LUT <- c("#938b72", "#E15759")

a<- ggplot(Cmb.avg, aes(x= Condition, y= Cmb.avg.lpd.db.mean.zs, color = Condition))
a+ geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.18, alpha=0.8, size=3, shape = 7)+
  scale_color_manual(values = M2.LUT)+
  theme_classic()+
  theme(text=element_text(size = 13, face = "plain"),
        axis.text = element_text(colour = "black"))+
  labs(title = "Mboat2 OE aging score" , x = "", y = "Mean Z Score", color = "")+
  facet_wrap(vars(Age), nrow = 1) +
  stat_pvalue_manual(stat.test, label = "p",
                     bracket.nudge.y = 0.1,
                     step.increase = 0.03, size = 3.5)+
  theme(legend.position= "none")

ggsave(filename = paste0("./Figure_Panels/Fig.4i.pdf"), width = 5, height = 5, useDingbats=FALSE)

med.sc <- Cmb.avg %>% 
  group_by(Age, Condition) %>% 
  summarise(med.score = median(Cmb.avg.lpd.db.mean.zs))
med.sc

aging.pval <- t.test(Cmb.avg$Cmb.avg.lpd.db.mean.zs[Cmb.avg$Age == "Young" & Cmb.avg$Condition == "Control"],Cmb.avg$Cmb.avg.lpd.db.mean.zs[Cmb.avg$Age == "Old" & Cmb.avg$Condition == "Control"], paired = F)
aging.pval$p.value
