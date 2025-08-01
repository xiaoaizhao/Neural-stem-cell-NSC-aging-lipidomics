## Mboat2 OE effect assessed by specific substrate and product level of target enzyme

rm(list = ls())
library(tidyverse)
library(ggpubr)
library(patchwork)
setwd(rstudioapi::getActiveProject())

M2.LUT <- c("#938b72", "#E15759")

load("./Output_Data/Mboat2_OE_LIPID.Rdata")

M2OE.l <- M2OE.Lipid %>% 
  rownames_to_column(var = "LipidIon") %>% 
  mutate(., LipidIon = ifelse(
    grepl("Ch\\+H\\-H2O|Cholesterol", LipidIon),
    "Cholesterol(0:0)+H-H20",
    LipidIon
  )) %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1))
## -------------------------------------------------------------------------------------------------------------------
## Mboat2 prefers aceyl-CoA 18:1 as substrate, but prefers LPE with these chains as acceptor in this order
## LPE 16:0 > 18:1 > alkenyl(crude) > 18:0
## -------------------------------------------------------------------------------------------------------------------

LPE_16_0 <- M2OE.l %>%
  filter(., Class == "LPE" & grepl("16:0", SideChain)) 

PE_16_0 <- M2OE.l %>%
  filter(., Class == "PE" & grepl("16:0", SideChain)) %>%
  group_by(., Sample) %>%
  summarise(., PE_16_0_all = sum(Conc_Int))

M.df <- left_join(PE_16_0, LPE_16_0, by = "Sample") %>%
  rename(., "LPE_16_0_conc_int" = "Conc_Int") %>%
  mutate(., PE_LPE_ratio = PE_16_0_all/LPE_16_0_conc_int) %>%
    mutate(Condition = case_when(
      grepl("_Mb2_OE", Sample) ~ "Mboat2 OE",
      grepl("_EGFP", Sample) ~ "Control"
    )) %>% 
    mutate(., Age = ifelse(grepl("^Y", Sample), "Young", "Old")) 

M.df$Age <- factor(M.df$Age, levels = c("Young", "Old"))
M.df$Condition <- factor(M.df$Condition, levels = c("Control", "Mboat2 OE"))


#Plot substrate level
a<- ggplot(M.df, aes(x= factor(Condition), y= LPE_16_0_conc_int, colour= Condition)) +
 geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.85,width = 0.2, size = 3)+
  theme_classic()+
  scale_color_manual(values = M2.LUT)+
  labs(title = "LPE(16:0)" , x = "", y = "Concentration (uM)", color = "")+
  theme(text=element_text(size = 12, face = "plain"),
        axis.text = element_text(colour = "black"))+
  stat_compare_means(aes(group = Condition), label = "p.format", paired = TRUE, method = "wilcox.test")+
  theme(legend.position= "none")+
  ylim(0, max(M.df$LPE_16_0_conc_int) + 0.0001) +
  theme(axis.text = element_text(colour = "black"))
ggsave(filename = "./Figure_Panels/EDFig.11i.left.pdf", width = 3.5, height = 5, useDingbats=FALSE)

#Plot product level
b<- ggplot(M.df, aes(x= factor(Condition), y= PE_16_0_all, colour= Condition)) +
 geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.85,width = 0.2, size = 3)+
  theme_classic()+
  scale_color_manual(values = M2.LUT)+
  labs(title = "PE(16:0)" , x = "", y = "Concentration (uM)", color = "")+
  theme(text=element_text(size = 12, face = "plain"))+
  stat_compare_means(aes(group = Condition), label = "p.format", paired = TRUE, method = "wilcox.test")+
  theme(legend.position= "none") +
  ylim(0, max(M.df$PE_16_0_all) + 0.01) +
  theme(axis.text = element_text(colour = "black"))
ggsave(filename = "./Figure_Panels/EDFig.11i.middle.pdf", width = 3.5, height = 5, useDingbats=FALSE)

#Plot ratio between substrate to product
c<- ggplot(M.df, aes(x= factor(Condition), y= PE_LPE_ratio, colour= Condition)) +
 geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.85,width = 0.2, size = 3)+
  theme_classic()+
  scale_color_manual(values = M2.LUT)+
  labs(title = "PE_LPE_ratio" , x = "", y = "Ratio - PE(16:0)/LPE(16:0)", color = "")+
  theme(text=element_text(size = 12, face = "plain"))+
  stat_compare_means(aes(group = Condition), label = "p.format", paired = TRUE, method = "wilcox.test")+
  theme(legend.position= "none") +
  theme(axis.text = element_text(colour = "black"))
ggsave(filename = "./Figure_Panels/EDFig.11i.right.pdf", width = 3.5, height = 5, useDingbats=FALSE)