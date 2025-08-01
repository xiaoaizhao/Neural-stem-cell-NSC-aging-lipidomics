## PCA analysis on lipidomics of qNSC with young plasma membrane lipid supplementation in vitro
setwd(rstudioapi::getActiveProject())
library(tidyverse)
rm(list=ls())
load("./Output_data/M2PM_Norm_Impt_log2_conc_275_noconc_48_lipids.Rdata")

Smp.Key <- read.csv("./Input_Data/March_Sample_list_071123_forR.csv", stringsAsFactors = F)
Smp.Key.e <- Smp.Key %>% 
  select(Sample.Name, Sample_ID)  %>% 
  rowwise() %>% 
  mutate(Exp = case_when(
    grepl("_EGFP|_Mb2_OE", Sample.Name) ~ "Mboat2_OE",
    grepl("_meth|lpd", Sample.Name) ~ "GPMV_Sup",
  )) %>% 
  filter(!Sample_ID == "XZ_45")

PM.sup.lg2 <- M2PM.Impt_norm_conc_no_conc_all %>% 
  select(-XZ_45) %>% 
  rename_at(vars(matches(Smp.Key.e$Sample_ID)), ~Smp.Key.e$Sample.Name) %>% 
  select(matches("_meth|lpd"))

my.lipids <-  as.data.frame(t(PM.sup.lg2))
df_pca <- prcomp(my.lipids, center = TRUE, scale. = T)

df_out <- as.data.frame(df_pca$x) %>% 
  rownames_to_column(., var = "Lpd_smple")

df_out.p <- df_out %>% 
  rowwise() %>% 
  mutate(GPMV_sup = case_when(
    grepl("_meth", Lpd_smple) ~ "Control",
    grepl("lpd", Lpd_smple) ~ "PM lipid supp"
  ))  %>%
  mutate(., Culture = substr(Lpd_smple, 1, str_locate(Lpd_smple, "_")-1)) %>% 
  mutate(., Age = ifelse(grepl("^Y", Lpd_smple), "Young", "Old")) %>% 
  relocate(c(Culture, Age, GPMV_sup, Lpd_smple), .after = Lpd_smple)

xplot <- "PC1"
yplot <- "PC2"

df_out.p <- df_out.p %>% 
  ungroup() 

df_out.p$Age <- factor(df_out.p$Age, levels = c("Young", "Old"))
df_out.p$GPMV_sup <- factor(df_out.p$GPMV_sup, levels = c("Control", "PM lipid supp"))
df_out.p <- as.data.frame(df_out.p)

pal4 <- c("cyan3", "magenta3")
p<-ggplot(df_out.p,aes(x=.data[[xplot]],y=.data[[yplot]]))

p+geom_point(size=3, alpha=0.8, stroke = 1, aes(color = Age, shape = GPMV_sup))+
  theme_classic()+
  xlab(paste0(xplot, ": ", format(summary(df_pca)$importance[2,xplot] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0(yplot, ": ", format(summary(df_pca)$importance[2,yplot] * 100,
                                  digits = 3), " % variance"))+
  scale_color_manual(values = pal4) +
  scale_shape_manual(values = c(5, 14)) +
  labs(title = "GPMV lipid Supplementation - rmv XZ_45")+
  theme(text=element_text(size = 13, face = "plain"))
 ggsave(filename = "./Figure_Panels/EDFig.13e.pdf", width = 5.5, height = 5, useDingbats=FALSE)
 
