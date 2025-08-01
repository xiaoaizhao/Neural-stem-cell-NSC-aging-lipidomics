# Batch 2 - GPMV supplementation lipidomic aging score plotting
# plot lipidomics aging score based on combined individual lipid signature and DB composition signature

## Calculate lipidomic aging score and compare across samples in Mboat2 OE experiments (from Batch #2)
## Aging lipidomic score is calculated as the mean z score across aging features 


##Result: 
#1. using the a more inclusive criteria (no filtering on number of studies detected, with p value <0.05),
#   Supplementation with young plasma membrane lipids did not lead to significant changes on lipidomic aging score in either young or old samples.

## --------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(ggpubr)
## Import lipidomic aging signature ------
load("./Output_Data/Meta_Lipid_signature.Rdata")
load("./Output_Data/Meta_DB_signature.Rdata")

## PM lipid supplementation samples
load("./Output_Data/PM_supp_DB_PCT.Rdata")
load("./Output_Data/PM_supp_LIPID.Rdata")

## Remove outlier XZ_45 - Y5_YlpdP5
PM.DB.df <- PM.sup.DB %>% 
  filter(!Sample == "XZ_45")#1633

PM.lipid.df <- PM.sup.Lipid %>% 
  select(-Y5_YlpdP5) #321 lipids and 23 samples

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
  filter(., Cla_DB %in% DB.hi.old) %>% 
  group_by(Sample.Name) %>%
  summarise(., mean_zs = mean(zscore, na.rm = TRUE)) %>%
  mutate(Condition = case_when(
    grepl("YlpdP5", Sample.Name) ~ "Young GPMV P5",
    grepl("YlpdP7", Sample.Name) ~ "Young GPMV P7",
    grepl("meth", Sample.Name) ~ "Control"
  ))%>%
  mutate(., Age = ifelse(grepl("^Y", Sample.Name), "Young", "Old")) %>% 
  rename(., "DB_mean_zs" = "mean_zs") 

## -------------------------------------------------------------------------------------------------------------------
lpd.zs <- PM.lipid.df %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample.Name", values_to = "Conc_Int") %>%
  group_by(LipidIon) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(Conc_Int))
  })

highold.mean.lpd.zs <- lpd.zs %>%
  mutate(., LipidID = substr(LipidIon, 1, str_locate(LipidIon, "\\)"))) %>% 
  filter(., LipidID %in% conc.lpd.hi.old) 

highold.mean.lpd.zs <- highold.mean.lpd.zs %>% 
  group_by(Sample.Name) %>%
  summarise(., mean_zs = mean(zscore, na.rm = TRUE)) %>%
  mutate(Condition = case_when(
    grepl("YlpdP5", Sample.Name) ~ "Young GPMV P5",
    grepl("YlpdP7", Sample.Name) ~ "Young GPMV P7",
    grepl("meth", Sample.Name) ~ "Control"
  ))%>%
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

Cmb.avg$Condition <- factor(Cmb.avg$Condition, levels = c("Control", "Young GPMV P5", "Young GPMV P7"))
Cmb.avg$Age <- factor(Cmb.avg$Age, levels = c("Young", "Old"))
## Normality test to determine if it's appropriate to use parametric statistical testing
## -------------------------------------------------------------------------------------------------------------------
## Test between control and treatment in each age group

norm.test <- list()
for (x in levels(Cmb.avg$Age)){
  
  d <- with(Cmb.avg, 
            Cmb.avg.lpd.db.mean.zs[Condition == "Control"] - Cmb.avg.lpd.db.mean.zs[Condition == "Young GPMV P7"])
  
  norm.test[[x]] <- shapiro.test(d)$p.value
}

norm.df <- bind_rows(norm.test)
print(norm.df)
## result: all p-values are above 0.05, does not reject the hypothesis of normality. Will perform t-test for 

## Pair-wise t-test for old and young samples between KO and control conditions
## -------------------------------------------------------------------------------------------------------------------
stat.test <- Cmb.avg %>%
  group_by(Age) %>%
  t_test(Cmb.avg.lpd.db.mean.zs ~ Condition, paired = FALSE, ref.group = "Control") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "Condition")
stat.test$p.adj <- format(stat.test$p.adj, digits = 2)

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
  theme(legend.position= "none")

ggsave(filename = "./Figure_Panels/EDFig.13d.pdf", width = 5, height = 5, useDingbats=FALSE)

