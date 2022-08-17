
##In vivo sorted lipidomics
##Number of lipid per class - pie chart
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

load("./Output_Data/Invivo_Norm_Impt_log2_all130_lipids.Rdata")
invivo.pie<- Impt_norm_conc_no_conc_all %>%
  mutate(., Class = substr(row.names(Impt_norm_conc_no_conc_all), 1, 
                           str_locate(row.names(Impt_norm_conc_no_conc_all), "\\(")-1)) %>%
  mutate(., lipid_cat = Class) %>%
  mutate(., lipid_cat = ifelse(grepl("Cer", Class), "Cer", lipid_cat)) %>%
  group_by(lipid_cat) %>%
  summarise(., nLipid = n())

load("./Input_Data/Pie_chart_color_LUT_full.name_2021-02-02.Rdata")

MG.LPI <- tibble(lipid_cat = c("MG", "LPI"),
                 Clr_list.1.18. = c("#9D7660", "#D7B5A6"),
                 Full_name = c("Monoacylglyceride (MG)", "Lysophosphatidylinositol(LPI)"))
lut.all <- bind_rows(a.full, MG.LPI)

invivo.df.clr <- left_join(invivo.pie, lut.all, by = "lipid_cat") %>%
  arrange(Full_name)

a <- ggplot(invivo.df.clr, aes(x="", y = nLipid, fill = factor(Full_name)))
a+  geom_bar(position="fill", stat="identity") + coord_polar("y")+
  theme_classic()+
  theme(axis.title=element_blank(), axis.line=element_blank(),
        axis.ticks=element_blank(), axis.text=element_blank(),
        plot.background = element_blank(), 
        plot.title=element_text(color="black",size=9,face="bold",hjust=0.5),
        strip.background = element_blank(), 
        strip.text.x = element_text(color = "transparent") )+
  scale_fill_manual(values = invivo.df.clr$Clr_list.1.18.) + 
  labs(title = "In vivo Number of lipids per class",fill = "")
ggsave(filename = paste0("./Figure_Panels/Fig_S2a.pdf"), width = 5, height = 5, useDingbats=FALSE)
