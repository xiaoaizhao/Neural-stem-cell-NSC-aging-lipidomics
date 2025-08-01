# Batch 2 - Mboat2 OE lipidomic aging score plotting
# plot lipidomics aging score based on combined individual lipid signature and DB composition signature

## Calculate lipidomic aging score and compare across samples in Mboat2 OE experiments (from Batch #2)
## Aging lipidomic score is calculated as the mean z score across aging features 


# Lipidomic aging score in all Mboat2 in vitro OE lipidomics

## --------------------------------------------------------------------
rm(list=ls())
library(rstatix)
library(tidyverse)
library(ggpubr)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")
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
  })

highold.mean.DB.zs <- DB.zs %>%
  filter(., Cla_DB %in% DB.hi.old) %>% 
  group_by(Sample.Name) %>%
  summarise(., mean_zs = mean(zscore, na.rm = TRUE)) %>%
  mutate(Condition = case_when(
    grepl("_Mb2_OE", Sample.Name) ~ "Mboat2 OE",
    grepl("_EGFP", Sample.Name) ~ "Control"
  )) %>% 
  mutate(., Age = ifelse(grepl("^Y", Sample.Name), "Young", "Old")) %>% 
  rename(., "DB_mean_zs" = "mean_zs") 

## -------------------------------------------------------------------------------------------------------------------
lpd.zs <- M2OE.Lipid %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample.Name", values_to = "Conc_Int") %>%
  group_by(LipidIon) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(Conc_Int))
  })

highold.mean.lpd.zs <- lpd.zs %>%
  mutate(., LipidID = substr(LipidIon, 1, str_locate(LipidIon, "\\)"))) %>% 
  filter(., LipidID %in% conc.lpd.hi.old) #4 out of 26 lipid features were detected in this dataset

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
      mutate(., Cmb.avg.lpd.db.mean.zs = mean(c(DB_mean_zs, Lpd_mean_zs)))
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
# stat.test$p.adj <- format(stat.test$p.adj, digits = 2)


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

ggsave(filename = "./Figure_Panels/Fig.4i.pdf", width = 5, height = 5, useDingbats=FALSE)

## -------------------------------------------------------------------------------------------------------------------
df.median <- Cmb.avg %>% 
  group_by(Age, Condition) %>% 
  summarise(Med.score = median(Cmb.avg.lpd.db.mean.zs))
df.median
