##Objective:
#1. Statistical testing between young and old samples within each KO condition
#2. Create a bar chart showing the class composition of lipids increase with age, and a separate chart for decrease with age.

##Note: Given the low sample number in this experiment (4 young and 4 old for each condition), there aren't many lipids pass stats with FDR <0.05, use Wilcoxon p value instead.
##For supplemental panel, plot lipid change with age with P-value <0.05
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_data/Exp2_Norm_Impt_backtoraw_all693_lipids.Rdata")

####Statistical testing between young and old samples within each KO condition####
all.lipid.df <- raw_int.exp2 %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Intensity") %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old")) %>%
  mutate(., SampleID = substr(Sample, 1, str_locate(Sample, "_")-1)) %>%
  mutate(., KO = substr(Sample, str_locate(Sample, "_"), nchar(Sample)))

all.lipid.df$KO <- factor(all.lipid.df$KO)

stat.all <- list()
for (ko.name in levels(all.lipid.df$KO)) {
  df.ko <- all.lipid.df %>%
    filter(., KO == ko.name)
  stat.all[[ko.name]] <- wilcox_stat(df.ko, Conc_Intensity, LipidIon)
  stat.all[[ko.name]] <- stat.all[[ko.name]] %>%
    mutate(., KO = ko.name)
}

Exp2.all.stat <- bind_rows(stat.all)
save(Exp2.all.stat, file = paste0("./Output_Data/Exp2_Lipid_stats_693_all_OvY.Rdata"))

###Bar chart for only significant lipids####============================================================#############################
###Bar chart for only significant lipids####============================================================#############################
####Create a bar chart showing the class composition of significant lipids increase with age, and a separate chart for decrease with age.####
load("./Output_Data/Exp2_Lipid_stats_693_all_OvY.Rdata")
stat_n <- Exp2.all.stat %>%
  filter(., KO == "_N") %>%
  rename(., "LipidIon" = "LipidIon...1") %>%
  select(., -LipidIon...4) %>%
  mutate(., FA = substr(LipidIon, str_locate(LipidIon, "\\("), str_locate(LipidIon, "\\)"))) %>%
  filter(., Wilcox_Pval < 0.05) %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., Class = ifelse(grepl("Cer", Class), "Cer", Class)) %>% 
  group_by(LipidIon) %>% 
  group_modify(~ {
    .x %>%
      mutate(., Linker = substr(FA, str_locate(FA, "e|p|o"), str_locate(FA, "e|p|o"))) %>% #Annotate ether lipid linker, to separate them from ester lipids
      mutate(., Linker = ifelse(!is.na(Linker), "ether lipids", ""))
  }) %>%
  ungroup()


up.pie <- stat_n %>%
  filter(., log2FC_OvY>0)  %>% # 74
  mutate(., Class = paste0(Class, 
                           "_", Linker)) %>% 
  ungroup() %>%
  mutate(., lipid_cat = substr(Class, 1, str_locate(Class, "\\_")-1)) %>%
  group_by(lipid_cat, Class) %>%
  summarise(., nLipid = n())


down.pie <- stat_n %>%  # 33
  filter(., log2FC_OvY<0)  %>% 
  mutate(., Class = paste0(Class, 
                           "_", Linker)) %>% 
  ungroup() %>%
  mutate(., lipid_cat = substr(Class, 1, str_locate(Class, "\\_")-1)) %>%
  group_by(lipid_cat, Class) %>%
  summarise(., nLipid = n())

##Load in color look up table to keep color schmem consistent throughout the paper for each lipid class
load("./Input_Data/Pie_chart_color_LUT_full.name_2021-02-02.Rdata")

####Up with age chart####
up.clr <- left_join(up.pie, a.full, by = "lipid_cat")

up.bar.pct <- up.clr %>%
  ungroup() %>% 
  mutate(., Up.all = sum(nLipid)) %>%
  mutate(., Class_Pct = nLipid/Up.all*100) %>%
  mutate(., Full_name.p = paste(Full_name, 
                                substr(Class, str_locate(Class, "_")+1, nchar(Class)),
                                sep = " ")) %>% 
  arrange(., Full_name.p)

up.bar.pct$Full_name.p <- factor(up.bar.pct$Full_name.p)


a <- ggplot(up.bar.pct, aes(x=fct_reorder(Full_name.p, Class_Pct), y=Class_Pct, 
                            fill=factor(Full_name.p)))
a+ geom_bar(width = 1, stat = "identity") +
  coord_flip()+
  scale_fill_manual(values = as.character(up.bar.pct$Clr_list.1.18.)) +
  theme_classic()+
  labs(y="% Total Number of Lipids", x="", title = "Increase With Age", fill = "")+
  theme(legend.position = "none")+
  theme(text=element_text(size = 13, face = "plain"),
  axis.text = element_text(colour = "black"))
ggsave(filename = paste0("./Figure_Panels/Fig_S1e_bar_UP.pdf"), width = 5, height = 5,
       useDingbats=FALSE)


####Down with age chart####
down.clr <- left_join(down.pie, a.full, by = "lipid_cat")

down.bar.pct <- down.clr %>%
  ungroup() %>% 
  mutate(., Down.all = sum(nLipid)) %>%
  mutate(., Class_Pct = nLipid/Down.all*100) %>%
  mutate(., Full_name.p = paste(Full_name, 
                                substr(Class, str_locate(Class, "_")+1, nchar(Class)),
                                sep = " ")) %>% 
  arrange(., Full_name.p)

down.bar.pct$Full_name.p <- factor(down.bar.pct$Full_name.p)


a <- ggplot(down.bar.pct, aes(x=fct_reorder(Full_name.p, Class_Pct), y=Class_Pct, 
                            fill=factor(Full_name.p)))
a+ geom_bar(width = 1, stat = "identity") +
  coord_flip()+
  scale_fill_manual(values = as.character(down.bar.pct$Clr_list.1.18.)) +
  theme_classic()+
  labs(y="% Total Number of Lipids", x="", title = "Decrease With Age", fill = "")+
  theme(legend.position = "none")+
  theme(text=element_text(size = 13, face = "plain"),
        axis.text = element_text(colour = "black"))
ggsave(filename = paste0("./Figure_Panels/Fig_S1e_bar_DOWN.pdf"), width = 5, height = 5,
       useDingbats=FALSE)
