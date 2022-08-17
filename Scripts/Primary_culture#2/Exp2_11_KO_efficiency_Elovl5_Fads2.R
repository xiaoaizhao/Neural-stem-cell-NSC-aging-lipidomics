## KO efficiency assessed by specific substrate and product level of target enzyme
rm(list=ls())
library(tidyverse)
library(ggpubr)

setwd(rstudioapi::getActiveProject())

load(file = "Output_Data/Exp2_FFA_concentration.rdata")
load(file = "./Output_Data/KO_LUT.Rdata")

## Elovl5 and Fads2 first, both target free fatty acid
## Elovl5 extends FFA 18:2 to 20:2, and 18:3 to 20:3
## -------------------------------------------------------------------------------------------------------------------
FFA.conc <- FFA.conc %>% 
  mutate(., `18:2/20:2 Ratio` = `18:2`/`20:2`) %>%
  mutate(., `18:3/20:3 Ratio` = `18:3`/`20:3`) %>%
  mutate(., `18:2/18:3 Ratio` = `18:2`/`18:3`) %>%
  mutate(., `20:2/20:3 Ratio` = `20:2`/`20:3`) 


KO.name <- "Elovl5"
ls <- c("Control", KO.name)
ctrl_elovl5_clr <- KO.LUT$Color[KO.LUT$KO %in% ls]

E.df <- FFA.conc %>%
  select(., c("Sample", "20:2", "20:3", "18:2", "18:3")) %>%
  pivot_longer(-Sample, names_to = "Fatty_acids", values_to = "Product_Conc") %>%
  mutate(., KO = substr(Sample, 3, 4)) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old")) %>%
  filter(., grepl("_N|_E", KO))

E.df$KO[grep("N", E.df$KO)] <- "Control"
E.df$KO[grep("E", E.df$KO)] <- "Elovl5 KO"

E.df$KO <- factor(E.df$KO, levels = c("Control", "Elovl5 KO"))

#combine all age, compare to ctrl

KO.classsum.delta<- E.df %>%
  mutate(., Fatty_acids1 = Fatty_acids) %>%
  group_by(., Fatty_acids1) %>%
  group_map(~ {
    ggplot(.x) + aes(KO, Product_Conc, color = KO)+
      geom_boxplot(outlier.shape = NA)+
      geom_jitter(alpha=0.7,width = 0.2, size = 3)+
      scale_colour_manual(values = ctrl_elovl5_clr)+
      theme_classic()+
      labs(title = paste0("Free fatty acid ", unique(.x$Fatty_acids)) , x = "", y = "Concentration (uM)", color = "")+
      theme(text=element_text(size = 12, face = "plain"),
            axis.text = element_text(colour = "black"))+
      stat_compare_means(aes(group = KO), label = "p.format", paired = TRUE)+
      theme(legend.position= "none")
      ggsave(filename = paste0("./Figure_Panels/Fig_S6b-c_", KO.name, "_",  
                               unique(.x$Fatty_acids), "_", ".pdf"), width = 3.5, height = 5, useDingbats=FALSE)
  })


## Fads2 add additional double bond. FFA 18:2 convert to 18:3, and FFA 20:2 convert to 20:3
## -------------------------------------------------------------------------------------------------------------------
KO.name <- "Fads2"
ls <- c("Control", KO.name)
ctrl_fads2_clr <- KO.LUT$Color[KO.LUT$KO %in% ls]

F.df <- FFA.conc %>%
  select(., c("Sample", "20:2", "20:3", "18:2", "18:3")) %>%
  pivot_longer(-Sample, names_to = "Fatty_acids", values_to = "Product_Conc") %>%
  mutate(., KO = substr(Sample, 3, 4)) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old")) %>%
  filter(., grepl("_N|_F", KO))

F.df$KO[grep("N", F.df$KO)] <- "Control"
F.df$KO[grep("F", F.df$KO)] <- "FADS2 KO"

F.df$KO <- factor(F.df$KO, levels = c("Control", "FADS2 KO"))

#combine all age, compare to ctrl
KO.classsum.delta<- F.df %>%
  mutate(., Fatty_acids1 = Fatty_acids) %>%
  group_by(., Fatty_acids1) %>%
  group_map(~ {
    ggplot(.x) + aes(KO, Product_Conc, color = KO)+
      geom_boxplot(outlier.shape = NA)+
      geom_jitter(alpha=0.7,width = 0.2, size = 3)+
      scale_colour_manual(values = ctrl_fads2_clr)+
      theme_classic()+
      labs(title = paste0("Free fatty acid ", unique(.x$Fatty_acids)) , x = "", y = "Concentration (uM)", color = "")+
      theme(text=element_text(size = 10, face = "bold"))+
      stat_compare_means(aes(group = KO), label = "p.format", paired = TRUE)+
      theme(legend.position= "none")
      ggsave(filename = paste0("./Figure_Panels/Fig_S6b-c_", KO.name, "_",  
                             unique(.x$Fatty_acids), "_", ".pdf"), width = 3.5, height = 5, useDingbats=FALSE)
  })

