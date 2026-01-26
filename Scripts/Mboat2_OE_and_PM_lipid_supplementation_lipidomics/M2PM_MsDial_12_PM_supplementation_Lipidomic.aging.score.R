### Lipidomic aging score in samples with plasma membrane lipid supplementation

## --------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(ggpubr)
## Import lipidomic aging signature ------
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/Meta_Lipid_signature.Rdata")
load("./Output_Data/Meta_DB_signature.Rdata")

## PM lipid supplementation samples
load("./Output_Data/PM_supp_DB_PCT.Rdata")
load("./Output_Data/PM_supp_LIPID.Rdata")

## Remove outlier XZ_45 - Y5_YlpdP5
PM.DB.df <- PM.sup.DB %>% 
  filter(!Sample == "Y5_YlpdP5")#1403

PM.lipid.df <- PM.sup.Lipid %>% 
  select(-Y5_YlpdP5) #397 lipids and 23 samples

## Get z score for each aging double bond composition features 
## -------------------------------------------------------------------------------------------------------------------
DB.zs <- PM.DB.df %>%
  mutate(., Cla_DB = paste0(Class, DB_num)) %>%
  group_by(Cla_DB) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(DB_Pct))
  })

highold.mean.DB.zs <- DB.zs %>%
  filter(., Cla_DB %in% conc.DB.hi.old) %>% 
  group_by(Sample) %>%
  summarise(., mean_zs = mean(zscore, na.rm = TRUE)) %>%
  mutate(Condition = case_when(
    grepl("YlpdP5", Sample) ~ "Young GPMV P5",
    grepl("YlpdP7", Sample) ~ "Young GPMV P7",
    grepl("meth", Sample) ~ "Control"
  ))%>%
  mutate(., Age = ifelse(grepl("^Y", Sample), "Young", "Old")) %>% 
  rename(., "DB_mean_zs" = "mean_zs") 

## ---------------------Lipid## -------------------------------------------------------------------------------------------------------------------
lpd.zs <- PM.lipid.df %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") %>%
  group_by(LipidIon) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(Conc_Int))
  }) %>% 
  ungroup()

highold.mean.lpd.zs <- MsD.lpd.rmv.abc(lpd.zs) %>%
  filter(., LipidIon %in% conc.lpd.hi.old) 

highold.mean.lpd.zs <- highold.mean.lpd.zs %>% 
  group_by(Sample) %>%
  summarise(., mean_zs = mean(zscore, na.rm = TRUE)) %>%
  mutate(Condition = case_when(
    grepl("YlpdP5", Sample) ~ "Young GPMV P5",
    grepl("YlpdP7", Sample) ~ "Young GPMV P7",
    grepl("meth", Sample) ~ "Control"
  ))%>%
  mutate(., Age = ifelse(grepl("^Y", Sample), "Young", "Old")) %>% 
  rename(., "Lpd_mean_zs" = "mean_zs") 

## -------------------------------------------------------------------------------------------------------------------
Cmb <- left_join(highold.mean.DB.zs, highold.mean.lpd.zs, by = c("Sample", "Condition", "Age")) 

Cmb.avg <- Cmb %>% 
  group_by(., Sample) %>% 
  group_modify(~{
    .x %>% 
      mutate(., Cmb.avg.lpd.db.mean.zs = mean(c(DB_mean_zs, Lpd_mean_zs)))
  }) 

Cmb.avg$Condition <- factor(Cmb.avg$Condition, levels = c("Control", "Young GPMV P5", "Young GPMV P7"))
Cmb.avg$Age <- factor(Cmb.avg$Age, levels = c("Young", "Old"))


## -------------------------------------------------------------------------------------------------------------------
### combine two supplemnetation experiments together
Cmb.avg.p <- Cmb.avg %>% 
  mutate(Cat = ifelse(Condition == "Control", "Control", "Young GPMV"))

Cmb.avg.p$Cat <- factor(Cmb.avg.p$Cat, levels = c("Control", "Young GPMV"))

pal4 <- c("cyan3", "magenta3")

a<- ggplot(Cmb.avg.p, aes(x= Cat, y= Cmb.avg.lpd.db.mean.zs))
a+ geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.18, alpha=0.8, stroke = 1, size = 3, aes(color = Age, shape = Cat))+
  scale_color_manual(values = pal4) +
  scale_shape_manual(values = c(5, 14)) +
  theme_classic()+
  theme(text=element_text(size = 13, face = "plain"),
        axis.text = element_text(colour = "black"))+
  labs(title = "PM lipid supp aging score" , x = "", y = "Mean Z Score", color = "")+
  facet_wrap(vars(Age), nrow = 1) +
  stat_compare_means(aes(group = Cat), label = "p.format", method = "wilcox.test") +
  theme(legend.position= "none") +
  ylim(min(Cmb.avg.p$Cmb.avg.lpd.db.mean.zs)-0.05, 1.3)

ggsave(filename = paste0("./Figure_Panels/EDFig.13d.pdf"), width = 5, height = 5, useDingbats=FALSE)

##==== statistical test between age in each condition====
a<- ggplot(Cmb.avg.p, aes(x= Age, y= Cmb.avg.lpd.db.mean.zs))
a+ geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.18, alpha=0.8, stroke = 1, size = 3, aes(color = Age, shape = Cat))+
  scale_color_manual(values = pal4) +
  scale_shape_manual(values = c(5, 14)) +
  theme_classic()+
  theme(text=element_text(size = 13, face = "plain"),
        axis.text = element_text(colour = "black"))+
  labs(title = "Compare age PM lipid supp aging score" , x = "", y = "Mean Z Score", color = "")+
  facet_wrap(vars(Cat), nrow = 1) +
  stat_compare_means(aes(group = Age), label = "p.format", method = "wilcox.test") +
  theme(legend.position= "none")

med.sc <- Cmb.avg.p %>% 
  group_by(Age, Cat) %>% 
  summarise(med.score = median(Cmb.avg.lpd.db.mean.zs))

med.sc
