
##Create Pie chart to show number of lipids per class in 2017 LC data
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

load(file = "./Output_data/Spike-in_norm_Mednorm_all_373_lipids_back_to_raw_int.Rdata")
lc.pie <- raw_int %>%
  mutate(., Class = substr(row.names(raw_int), 1, str_locate(row.names(raw_int), "\\(")-1)) %>%
  mutate(., lipid_cat = Class) %>%
  mutate(., lipid_cat = ifelse(grepl("Cer", Class), "Cer", lipid_cat)) %>% ##combine all ceramide classes together
  group_by(lipid_cat) %>%
  summarise(., nLipid = n())

##apply the same color LUT for plotting
load("./Input_Data/Pie_chart_color_LUT_full.name_2021-02-02.Rdata")

##manually add MG and LPI since these two classes were not detected in other datasets, therefore not in the color LUT
MG.LPI <- tibble(lipid_cat = c("MG", "LPI"),
                 Clr_list.1.18. = c("#9D7660", "#D7B5A6"),
                 Full_name = c("Monoacylglyceride (MG)", "Lysophosphatidylinositol(LPI)"))
lut.all <- bind_rows(a.full, MG.LPI)

df.clr <- left_join(lc.pie, lut.all, by = "lipid_cat") %>%
  arrange(Full_name)

a <- ggplot(df.clr, aes(x="", y = nLipid, fill = factor(Full_name)))
a+  geom_bar(position="fill", stat="identity") + coord_polar("y")+
  theme_classic()+
  theme(axis.title=element_blank(), axis.line=element_blank(),
        axis.ticks=element_blank(), axis.text=element_blank(),
        plot.background = element_blank(), 
        plot.title=element_text(color="black",size=12,face="plain",hjust=0.5),
        strip.background = element_blank(), 
        strip.text.x = element_text(color = "transparent") )+
  scale_fill_manual(values = df.clr$Clr_list.1.18.) + 
  labs(title = "Number of lipids per class",fill = "")
ggsave(filename = paste0("./Figure_Panels/Fig_S1a.pdf"), width = 5, height = 5, useDingbats=FALSE)
