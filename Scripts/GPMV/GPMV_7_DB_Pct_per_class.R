
##Double bond relative abundance per class
##Steps:
##1. Calculate class sum
##2. Take accumulative intensity of each double bond for each class divided by class sum intensity
##3. Plot double bond abundance difference between age groups only in classes that have at least 10 lipids identified
##4. Plot Phospholipid order index based on sum concentration of (PE+SM)/PC.
rm(list=ls())
library(tidyverse)
library(ggpubr)
setwd(rstudioapi::getActiveProject())

##Get class sum#############################################################################################################
load("./Output_Data/GPMV_Norm_Impt_backtoraw_all565_lipids.Rdata")
all.lipid <- raw_int.GPMV %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") %>% #4476/373=12 samples
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1))

GPMV.classsum <- all.lipid %>%
  group_by(Sample, Class) %>%
  summarise(., Class_sum = sum(Conc_Int))

save(GPMV.classsum, file = "./Output_Data/GPMV_ClassSum_565_lipids.Rdata")

##Get composition PCT for DB in each class################################################################################
load(file = "./Output_data/GPMV_DB_by_Class.Rdata")
DB_by_class_GPMV$Class <- as.character(DB_by_class_GPMV$Class)
GPMV.classsum$Class <- as.character(GPMV.classsum$Class)

GPMV_DB <- left_join(DB_by_class_GPMV, GPMV.classsum) %>%
  mutate(., DB_Pct = Sum_DB/Class_sum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(GPMV_DB, file = paste0("./Output_Data/GPMV_DB_PCT_all_samples.Rdata"))

##DB pct age comparison ###################################################################################################
OY_mean <- GPMV_DB %>%
  group_by(Class, DB_num, Age) %>%
  summarise(., DB_Age_avg = mean(DB_Pct)) 

OY_FC <- OY_mean %>%
  group_by(Class, DB_num) %>%
  group_modify(~ {
    .x %>%
      mutate(., log2FC_OvY = log2(DB_Age_avg/DB_Age_avg[Age == "Young"]))
  }) %>%
  filter(., log2FC_OvY !=0) %>%
  select(., -matches("Age|Age_Avg"))

####Filter to only plot class that contains at least 10 lipids###########################################################
Lipid_numLC <- raw_int.GPMV %>%
  rownames_to_column(., var = "LipidIon") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  group_by(Class) %>%
  summarise(., nLipid = n()) %>%
  filter(.,nLipid>=10 )

GPMV_DB_filter <- OY_FC %>%
  filter(., Class %in% Lipid_numLC$Class)

##plot with class list which has at least 10 lipids detected in this experiments#########################################
##plot with class list which has at least 10 lipids detected in this experiments#########################################
p <- ggplot(data =  GPMV_DB_filter, aes(x= Class, y= DB_num, fill= log2FC_OvY)) 
p +  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0)+
  labs(fill= "Log2 FC")+
  labs(title = "GPMV Old vs. Young DB abundance", x = "Class", y = "Double Bond Number", color = "")+
  theme(text=element_text(size = 13, face = "plain"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(paste0("./Figure_Panels/Fig_3d.pdf"), width = 5, height = 5, useDingbats=FALSE)


##Plot Phospholipid Order Index based on class sum of (PE+SM)/PC#########################################################
##Plot Phospholipid Order Index based on class sum of (PE+SM)/PC#########################################################
load("./Output_Data/GPMV_ClassSum_565_lipids.Rdata")
sum.df <- GPMV.classsum %>%
  pivot_wider(., names_from = "Class", values_from = "Class_sum") %>%
  mutate(., POI = (PE+SM)/PC) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

sum.df$Age <- factor(sum.df$Age, levels = c("Young", "Old"))
c <- ggplot(sum.df, aes(Age, POI, color = Age))
c+  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.1, size = 3, alpha = 0.7)+
  theme_classic()+
  stat_compare_means(aes(group = Age), label = "p.format")+
  labs(title = "GPMV Phospho Order Index", y = "(PE+SM)/PC")+
  scale_color_manual(values = c("darkgoldenrod", "maroon"))+
  theme(legend.position= "none")+
  theme(text=element_text(size = 13, face = "plain"))
ggsave(paste0("./Figure_Panels/Fig_3f.pdf"), width = 3.5, height = 5, useDingbats=FALSE)

sum.df.med <- sum.df %>% 
  group_by(Age) %>% 
  summarise(., MedianPOI = median(POI))

sum.df.med
