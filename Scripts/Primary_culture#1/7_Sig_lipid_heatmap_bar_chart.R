##Objective:
#1. Statistical testing between young and old samples of qNSC
#2. Create heatmap on significant lipid with age on 2017 LC data
#3. Creat a bar chart showing the class composition of significant lipids increase with age, and a separate chart for decrease with age.
rm(list=ls())
library(tidyverse)
library(ComplexHeatmap)
library(GetoptLong)
library(stringi)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")
load(file = "./Output_data/Spike-in_norm_Mednorm_all_373_lipids_back_to_raw_int.Rdata")

##Subset quiescent NSC samples for statistical testing
qui.lipid <- raw_int %>%
  rownames_to_column(., var = "LipidIon")  %>%
  select(., matches("Qui|LipidIon")) %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Intensity") #4476/373=12 samples

all.lipid.df <- qui.lipid %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))


###Statistical test using wilcoxon rank sum test####
stat.all <- wilcox_stat(all.lipid.df, Intensity, LipidIon)
Qui_only_2017.stat <- stat.all
save(Qui_only_2017.stat, file = paste0("./Output_data/Qui_sample_OvY_stats_2017_LC.Rdata"))

###########
load("./Output_Data/Qui_sample_OvY_stats_2017_LC.Rdata")
load(file = "./Output_data/Spike-in_norm_Mednorm_all_373_lipids_back_to_raw_int.Rdata")
q.sig <- Qui_only_2017.stat %>%
  filter(., Padj<0.05) #95 lipids
sig.list <- q.sig$LipidIon...1

#Change side separation from "/" to "_", and remove ion info
#Annotate ether lipid linker. This will be used to plot ether lipids separately in the bar graph.
sig.heatmap.df <- raw_int %>%
  select(., matches("Quiescent")) %>%
  rownames_to_column(., var = "LipidIon") %>%
  filter(., LipidIon %in% sig.list) %>%
  mutate(., LipidIon = str_replace_all(LipidIon, "/", "_")) %>%
  mutate(., LipidIon = substr(LipidIon, 1, str_locate(LipidIon, "\\)"))) %>%
  mutate(., FA = substr(LipidIon, str_locate(LipidIon, "\\("), str_locate(LipidIon, "\\)"))) %>%
  mutate(., LipidIon1 = LipidIon) %>%
  group_by(LipidIon) %>%
  group_modify(~ {
    .x %>%
    mutate(., Linker = substr(FA, str_locate(FA, "d|e|p|o|t"), str_locate(FA, "d|e|p|o|t"))) %>%
    mutate(., Linker = ifelse(is.na(Linker), "", Linker)) %>%
    mutate(., LipidIon.P = paste0(substr(LipidIon1, 1, str_locate(LipidIon1, "\\(")-1), Linker, FA))
  }) %>%
  ungroup() 

###Heatmap with complex heatmap package####============================================================#############################
#Generate z score matrix first
z.sig.df <- sig.heatmap.df %>%
  select(., -c("FA", "Linker", "LipidIon", "LipidIon1")) %>%
  rename(., "LipidIon" = "LipidIon.P") %>%
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Intensity") %>%
  group_by(LipidIon) %>%
  group_modify(~{
    .x %>%
      mutate(., zscore = scale(Intensity))
  })

z.sig.plot <- z.sig.df %>%
  select(., -Intensity) %>%
  pivot_wider(names_from = Sample, values_from = zscore) %>%
  column_to_rownames(., var = "LipidIon")

#Get lipid names from 2 largest clusters (e.g. one cluster goes up with age, the other goes down with age)
aa <- Heatmap(z.sig.plot, km=2)
aa <- draw(aa)

clus.order <- row_order(aa)

clus1 <- z.sig.plot %>%
  rownames_to_column(., var = "LipidIon") %>%
  rowid_to_column(., var = "rowID") %>%
  column_to_rownames(., var = "LipidIon") %>%
  filter(., rowID %in% clus.order$`1`) %>%
  select(., -rowID)

clus2 <- z.sig.plot %>%
  rownames_to_column(., var = "LipidIon") %>%
  rowid_to_column(., var = "rowID") %>%
  column_to_rownames(., var = "LipidIon") %>%
  filter(., rowID %in% clus.order$`2`) %>%
  select(., -rowID)

clus.ord.mtx <- bind_rows(clus2, clus1)

pdf(qq("./Figure_Panels//Fig_1c.pdf"), width =4, height =8)
HM.2 <- Heatmap(clus.ord.mtx, name = "mat", 
                cluster_rows = FALSE,
                column_title = "Lipids with significant change with age in qNSCs",
                row_names_gp = gpar(fontsize = 6.2),
                )
draw(HM.2)
dev.off()

###Pie chart for only significant lipids####============================================================#############################
###Pie chart for only significant lipids####============================================================#############################
load("./Output_Data/Qui_sample_OvY_stats_2017_LC.Rdata")
pie <- Qui_only_2017.stat %>%
  rename(., "LipidIon" = "LipidIon...1") %>%
  select(., -LipidIon...4) %>%
  mutate(., FA = substr(LipidIon, str_locate(LipidIon, "\\("), str_locate(LipidIon, "\\)"))) %>%
  mutate(., LipidIon1 = LipidIon) %>%
  group_by(LipidIon) %>%
  group_modify(~ {
    .x %>%
      mutate(., Linker = substr(FA, str_locate(FA, "e|p|o"), str_locate(FA, "e|p|o"))) %>% #Annotate ether lipid linker, to separate them from ester lipids
      mutate(., Linker = ifelse(!is.na(Linker), "ether lipids", ""))
  }) %>%
  ungroup()

up.pie <- pie %>%
  filter(., Padj<0.05) %>%
  filter(., log2FC_OvY>0)  %>% # 40
  mutate(., Class = paste0(substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1), 
                           "_", Linker)) %>% 
  ungroup() %>%
  mutate(., lipid_cat = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  group_by(lipid_cat, Class) %>%
  summarise(., nLipid = n())

down.pie <- pie %>%
  filter(., Padj<0.05) %>%
  filter(., log2FC_OvY<0)  %>% #55
  mutate(., Class = paste0(substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1),
                           "_", Linker)) %>% #Annotate ether lipid linker, to separate them from ester lipids
  mutate(., Class = ifelse(grepl("TG", Class), "TG_", Class)) %>%
  ungroup() %>%
  mutate(., lipid_cat = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>% 
  group_by(lipid_cat, Class) %>%
  summarise(., nLipid = n()) 
  

##Load in color look up table to keep color schmem consistent throughout the paper for each lipid class
load("./Input_Data/Pie_chart_color_LUT_full.name_2021-02-02.Rdata")
up.clr <- left_join(up.pie, a.full, by = "lipid_cat")

####Up with age chart####
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
  theme(text=element_text(size = 13, face = "plain"))
ggsave(filename = paste0("./Figure_Panels/Fig_1d_bar_UP.pdf"), width = 7, height = 5, 
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

bar <- ggplot(down.bar.pct, aes(x=fct_reorder(Full_name.p, Class_Pct), y=Class_Pct, fill=factor(Full_name.p)))
bar + geom_bar(width = 1, stat = "identity") +
  coord_flip()+
  scale_fill_manual(values = as.character(down.bar.pct$Clr_list.1.18.))+
  theme_classic()+
  labs(y="% Total Number of Lipids", x="", title = "Decrease With Age", fill = "")+
  theme(legend.position = "none")+
  theme(text=element_text(size = 13, face = "plain"))
ggsave(filename = paste0("./Figure_Panels/Fig_1d_bar_DOWN.pdf"), width = 7, height = 5, useDingbats=FALSE)
