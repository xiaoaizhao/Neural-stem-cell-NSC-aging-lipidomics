#### Plot DB aging signature side-by-side with lipid aging signature
rm(list=ls())
library(tidyverse)
library(patchwork)
setwd(rstudioapi::getActiveProject())

##### Double bond signature first
load("./Output_Data/Meta_DoubleBond_4_studies.Rdata")
load("./Output_Data/Meta_DB_signature_hi_in_Old_4_studies.Rdata")
load("./Output_Data/Meta_DB_signature_hi_in_Yng_4_studies.Rdata")

summaryCastDT.2 <- summaryCastDT.DB %>%
  mutate(., dataset = "Summary") %>%
  rename(., 'es_g' = 'summary') %>%
  rename(., 'se_g' = 'se_summary')

Old.DB.Sig <- summaryCastDT.2 %>% 
  filter(., Cla_DB %in% DB.hi.old)

Old.DB.Sig$Cla_DB <- fct_reorder(Old.DB.Sig$Cla_DB, Old.DB.Sig$es_g)

DB.frst <- ggplot(Old.DB.Sig, aes(x = es_g, y = Cla_DB, xmin = es_g - se_g, xmax = es_g + se_g))+
  geom_point(size = 3, colour = "maroon") +
  geom_errorbarh(height=.2, colour = "maroon") +
  theme_classic()+
  theme(text=element_text(size = 12, face = "plain", colour = "black"),
        axis.text = element_text(colour = "black"))+
  labs(title = "Meta analysis on Double Bond Composition" , x = "Effect Size", y = "")+
  theme(legend.position= "none")+
  geom_vline(xintercept = 0, linetype = "dashed")


##### Lipid signature 
load("./Output_Data/Meta_lipids_4_studies.Rdata")
load("./Output_Data/Meta_Lipid_signature_hi_in_Old_4_studies.Rdata")
load("./Output_Data/Meta_Lipid_signature_hi_in_Yng_4_studies.Rdata")

summaryCastDT.2 <- summaryCastDT %>%
  mutate(., dataset = "Summary") %>%
  rename(., 'es_g' = 'summary') %>%
  rename(., 'se_g' = 'se_summary')

Old.Lpd.Sig <- summaryCastDT.2 %>% 
  filter(., LipidIon %in% lpd.hi.old) %>% 
  mutate(., LipidNoIon = substr(LipidIon, 1, str_locate(LipidIon, "\\)")))

Old.Lpd.Sig$LipidNoIon <- fct_reorder(Old.Lpd.Sig$LipidNoIon, Old.Lpd.Sig$es_g)

lpd.frst <- ggplot(Old.Lpd.Sig, aes(x = es_g, y = LipidNoIon, xmin = es_g - se_g, xmax = es_g + se_g)) +
  geom_point(size = 3, colour = "maroon") +
  geom_errorbarh(height=.2, colour = "maroon") +
  theme_classic()+
  theme(text=element_text(size = 12, face = "plain", colour = "black"),
        axis.text = element_text(colour = "black"))+
  labs(title = "Meta analysis on individual lipids" , x = "Effect Size", y = "")+
  theme(legend.position= "none")+
  geom_vline(xintercept = 0, linetype = "dashed")

##### Combine two plots together
DB.frst + lpd.frst
ggsave(filename = "./Figure_Panels/Fig_3e.pdf", width = 5, height = 5, useDingbats=FALSE)
