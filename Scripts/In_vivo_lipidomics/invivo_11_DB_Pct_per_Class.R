
##Double bond relative abundance per class
##Steps:
##1. Calculate class sum
##2. Take accumulative intensity of each double bond for each class divided by class sum intensity
##The above 2 steps are identical to Primary culture LC-MS data
##3. Plot double bond abundance difference between age groups in all classes rather than filter by class size - 
##This is because in vivo data has significantly less number of detected lipids. With a class size filter >10 there will only be 4 classes left.
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Pre-processing_functions.R")

##Get class sum####
load("./Output_Data/Invivo_Norm_Impt_backtoraw_all130_lipids.Rdata")
all.lipid <- raw_int.invivo %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") %>% #4476/373=12 samples
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1))

Invivo.classsum <- all.lipid %>%
  group_by(Sample, Class) %>%
  summarise(., Class_sum = sum(Conc_Int))

save(Invivo.classsum, file = "./Output_Data/Invivo_Class.sum_ConcInt.Rdata")
##Get composition PCT for DB in each class####
load(file = "./Output_data/Invivo_DB_by_Class_qNSC.Rdata")
DB_by_class_invivo$Class <- as.character(DB_by_class_invivo$Class)
Invivo.classsum$Class <- as.character(Invivo.classsum$Class)

Invivo.DB.PCT <- left_join(DB_by_class_invivo, Invivo.classsum) %>%
  mutate(., DB_Pct = Sum_DB/Class_sum) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

save(Invivo.DB.PCT, file = "./Output_Data/Invivo__DB_PCT_by_Class.Rdata")

##DB pct age comparison

invivo.OY_FC <- DB.FC(Invivo.DB.PCT)

##plot data from all classes #########################################
##plot data from all classes #########################################
p <- ggplot(data =  invivo.OY_FC, aes(x= Class, y= DB_num, fill= log2FC_OvY)) 
p +  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0)+
  labs(fill= "Log2 FC")+
  labs(title = "In vivo qNSC Old vs. Young DB abundance", x = "Class", y = "Double Bond Number", color = "")+
  theme(text=element_text(size = 13, face = "plain"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(paste0("./Figure_Panels/Fig_S2c.pdf"), width = 5, height = 5, useDingbats=FALSE)

