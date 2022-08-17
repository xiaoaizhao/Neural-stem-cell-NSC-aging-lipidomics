
##Double bond relative abundance per class
##Steps:
##1. Calculate class sum
##2. Take accumulative intensity of each double bond for each class divided by class sum intensity
##3. Plot double bond abundance difference between age groups only in classes that have at least 10 lipids identified
rm(list=ls())
library(tidyverse)
source("./Scripts/Function_scripts/Effect_size_functions.R")
setwd(rstudioapi::getActiveProject())

##Get class sum#############################################################################################################
load("./Output_Data/Exp2_Norm_Impt_backtoraw_all693_lipids.Rdata")
all.lipid <- raw_int.exp2 %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Intensity") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., KO = substr(Sample, str_locate(Sample, "_"), nchar(Sample))) %>%
  mutate(., ID = substr(Sample, 1, str_locate(Sample, "_")-1))

Exp2.classsum <- all.lipid %>%
  group_by(KO, ID, Class) %>%
  summarise(., Class_sum = sum(Conc_Intensity))
save(Exp2.classsum, file = paste0("./Output_Data/Exp2_class_sum_693_lipids.Rdata"))

##Get composition PCT for DB in each class################################################################################
load("./Output_Data/Exp2_DB_by_Class.Rdata")
Exp2_DB_by_class$Class <- as.character(Exp2_DB_by_class$Class)
Exp2.classsum$Class <- as.character(Exp2.classsum$Class)

Exp2.classsum <- Exp2.classsum %>%
  mutate(., Sample = paste0(ID, KO))

Exp2_DB <- left_join(Exp2_DB_by_class, Exp2.classsum) %>%
  mutate(., DB_Pct = Sum_DB/Class_sum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(Exp2_DB, file = paste0("./Output_Data/Exp2_DB_PCT_all_samples.Rdata"))

##DB pct age comparison in each KO group##################################################################################
Exp2_DB$KO <- factor(Exp2_DB$KO)
Exp2_DB <- Exp2_DB %>%
  mutate(., Cla_DB = paste0(Class, DB_num))

DB.fc.all <- list()

for (ko.name in levels(Exp2_DB$KO)) {
  df.ko <- Exp2_DB %>%
    filter(., KO == ko.name)
  DB.fc.all[[ko.name]] <- wilcox_stat(df.ko, DB_Pct, Cla_DB)
  DB.fc.all[[ko.name]] <- DB.fc.all[[ko.name]] %>%
    mutate(., KO = ko.name)
}

Exp2_DBpct_all_clas_w_stats <- bind_rows(DB.fc.all)

####Filter to only plot class that contains at least 10 lipids###########################################################
load("./Output_Data/Exp2_Norm_Impt_backtoraw_all693_lipids.Rdata")
Lipid_numLC <- raw_int.exp2 %>%
  rownames_to_column(., var = "LipidIon") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  group_by(Class) %>%
  summarise(., nLipid = n()) %>%
  filter(.,nLipid>=10 )

Exp2_DB_filter <- Exp2_DBpct_all_clas_w_stats %>%
  mutate(., Class = substr(Cla_DB...1, 1, str_locate(Cla_DB...1, ":")-1)) %>%
  mutate(., DB_num = substr(Cla_DB...1, str_locate(Cla_DB...1, ":"), nchar(Cla_DB...1))) %>%
  filter(., Class %in% Lipid_numLC$Class)

##plot with class list which has at least 10 lipids detected in this experiments#########################################
##plot with class list which has at least 10 lipids detected in this experiments#########################################
##Plot Non-targeting control sample only
df.plot <- Exp2_DB_filter %>%
  filter(., KO == "_N")

a <- ggplot(df.plot, aes(x= Class, y= DB_num, fill= log2FC_OvY))
a+geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0)+
  labs(fill= "Log2 FC")+
  labs(title = paste0("Non-targeting Old vs. Young DB abundance"), x = "Class", y = "Double Bond Number", color = "")+
  theme(text=element_text(size = 13, face = "plain"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(filename = "./Figure_Panels/Fig_S1i.pdf", width = 5, height = 5, useDingbats=FALSE)
