## KO efficiency assessed by specific substrate and product level of target enzyme
rm(list=ls())
library(tidyverse)
library(ggpubr)

setwd(rstudioapi::getActiveProject())

load(file = "./Output_Data/KO_LUT.Rdata")
load(file = "./Output_Data/Exp2_Norm_Impt_backtoraw_all693_lipids.Rdata")

KO.df <- raw_int.exp2 %>%
  rownames_to_column(., var = "LipidIon") %>% 
  pivot_longer(-LipidIon, names_to = "Sample", values_to = "Conc_Int") %>%
  mutate(., Class = substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1)) %>%
  mutate(., SideChain = substr(LipidIon, str_locate(LipidIon, "\\(")+1, str_locate(LipidIon, "\\)")-1))

## -------------------------------------------------------------------------------------------------------------------
## Mboat2 prefers aceyl-CoA 18:1 as substrate, but prefers LPE with these chains as acceptor in this order
## LPE 16:0 > 18:1 > alkenyl(crude) > 18:0
## -------------------------------------------------------------------------------------------------------------------
KO.name <- "Mboat2"
ls <- c("Control", KO.name)
ctrl_mboat2_clr <- KO.LUT$Color[KO.LUT$KO %in% ls]


LPE_16_0 <- KO.df %>%
  filter(., Class == "LPE" & grepl("16:0", SideChain)) 

PE_16_0 <- KO.df %>%
  filter(., Class == "PE" & grepl("16:0", SideChain)) %>%
  mutate(., KO = substr(Sample, 3, 4)) %>%
  group_by(., Sample) %>%
  summarise(., PE_16_0_all = sum(Conc_Int))

M.df <- left_join(PE_16_0, LPE_16_0, by = "Sample") %>%
  rename(., "LPE_16_0_conc_int" = "Conc_Int") %>%
  mutate(., LPE_PE_ratio = LPE_16_0_conc_int/PE_16_0_all) %>%
  mutate(., KO = substr(Sample, 3, 4)) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old")) %>%
  filter(., grepl("_N|_M", KO))

M.df$KO[grep("N", M.df$KO)] <- "Control"
M.df$KO[grep("M", M.df$KO)] <- "Mboat2 KO"

M.df$Age <- factor(M.df$Age, levels = c("Young", "Old"))
M.df$KO <- factor(M.df$KO, levels = c("Control", "Mboat2 KO"))


#Plot substrate level
a<- ggplot(M.df, aes(x= factor(KO), y= LPE_16_0_conc_int, colour= KO))
a+ geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.7,width = 0.2, size = 3)+
  theme_classic()+
  scale_colour_manual(values = ctrl_mboat2_clr)+
  labs(title = "LPE(16:0)" , x = "", y = "Concentration (uM)", color = "")+
  theme(text=element_text(size = 12, face = "plain"),
        axis.text = element_text(colour = "black"))+
  stat_compare_means(aes(group = KO), label = "p.format", paired = TRUE)+
  theme(legend.position= "none")
ggsave(filename = paste0("./Figure_Panels/Fig_S6b-c_", KO.name, "_",  
                         "substrate", ".pdf"), width = 3.5, height = 5, useDingbats=FALSE)

#Plot product level
a<- ggplot(M.df, aes(x= factor(KO), y= PE_16_0_all, colour= KO))
a+ geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.7,width = 0.2, size = 3)+
  theme_classic()+
  scale_colour_manual(values = ctrl_mboat2_clr)+
  labs(title = "PE(16:0)" , x = "", y = "Concentration (uM)", color = "")+
  theme(text=element_text(size = 10, face = "bold"))+
  stat_compare_means(aes(group = KO), label = "p.format", paired = TRUE)+
  theme(legend.position= "none")
ggsave(filename = paste0("./Figure_Panels/Fig_S6b-c_", KO.name, "_",  
                         "product", ".pdf"), width = 3.5, height = 5, useDingbats=FALSE)

## -------------------------------------------------------------------------------------------------------------------
## Pla2g4e seems have very weak phospholipase (PLA1 and PLA2 activity). It also have acyltransferase activity (generate NAPE)
## NAPE: use PC (mostly) and PE as O-acyl donor and PE as acyl acceptor to generate NAPE.
## PC/ LPC 16:0, this is the reaction catalized by Pla2g4e. Use this ratio to assess Pla2g4e KO efficiency
## -------------------------------------------------------------------------------------------------------------------
KO.name <- "Pla2g4e"
ls <- c("Control", KO.name)
ctrl_pla_clr <- KO.LUT$Color[KO.LUT$KO %in% ls]


LPC_16_0 <- KO.df %>%
  filter(., Class == "LPC" & grepl("16:0", SideChain)) %>%
  mutate(., KO = substr(Sample, 3, 4)) %>%
  group_by(., Sample) %>%
  summarise(., LPC_16_0_all = sum(Conc_Int))

PC_16_0 <- KO.df %>%
  filter(., Class == "PC" & grepl("16:0", SideChain)) %>%
  mutate(., KO = substr(Sample, 3, 4)) %>%
  group_by(., Sample) %>%
  summarise(., PC_16_0_all = sum(Conc_Int))

P.df <- left_join(PC_16_0, LPC_16_0, by = "Sample") %>%
  mutate(., PC_LPC_ratio = PC_16_0_all/LPC_16_0_all) %>%
  mutate(., KO = substr(Sample, 3, 4)) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old")) %>%
  filter(., grepl("_N|_P", KO))

P.df$KO[grep("N", P.df$KO)] <- "Control"
P.df$KO[grep("P", P.df$KO)] <- "Pla2g4e KO"


P.df$Age <- factor(P.df$Age, levels = c("Young", "Old"))
P.df$KO <- factor(P.df$KO, levels = c("Control", "Pla2g4e KO"))

#Plot substrate level
a<- ggplot(P.df, aes(x= factor(KO), y= PC_16_0_all, colour= KO))
a+ geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.7,width = 0.2, size = 3)+
  theme_classic()+
  scale_colour_manual(values = ctrl_pla_clr)+
  labs(title = "PC(16:0)" , x = "", y = "Concentration (uM)", color = "")+
  theme(text=element_text(size = 12, face = "plain"),
        axis.text = element_text(colour = "black"))+
  stat_compare_means(aes(group = KO), label = "p.format", paired = TRUE)+
  theme(legend.position= "none")
ggsave(filename = paste0("./Figure_Panels/Fig_S6b-c_", KO.name, "_",  
                         "substrate", ".pdf"), width = 3.5, height = 5, useDingbats=FALSE)

#Plot product level
a<- ggplot(P.df, aes(x= factor(KO), y= LPC_16_0_all, colour= KO))
a+ geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.7,width = 0.2, size = 3)+
  theme_classic()+
  scale_colour_manual(values = ctrl_pla_clr)+
  labs(title = "LPC(16:0)" , x = "", y = "Concentration (uM)", color = "")+
  theme(text=element_text(size = 10, face = "bold"))+
  stat_compare_means(aes(group = KO), label = "p.format", paired = TRUE)+
  theme(legend.position= "none")
ggsave(filename = paste0("./Figure_Panels/Fig_S6b-c_", KO.name, "_",  
                         "product", ".pdf"), width = 3.5, height = 5, useDingbats=FALSE)
