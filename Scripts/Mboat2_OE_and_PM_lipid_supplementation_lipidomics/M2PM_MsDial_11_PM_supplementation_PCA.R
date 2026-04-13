### PCA analysis on lipidomics of qNSC with young plasma membrane lipid supplementation in vitro
setwd(rstudioapi::getActiveProject())
library(tidyverse)
rm(list=ls())

load("./Output_data/M2PM_MsD.Norm_Impt_log2_conc_397_lipids.Rdata")

PM.sup.Lipid <- M2PM.log2.impt.norm.conc.MsD %>% 
  select(matches("_meth|lpd")) # 24 samples

## Remove outlier XZ_45 - Y5_YlpdP5
PM.sup.Lipid <- PM.sup.Lipid %>% 
  select(-Y5_YlpdP5)

my.lipids <-  as.data.frame(t(PM.sup.Lipid))
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
  labs(title = "GPMV lipid Supplementation")+
  theme(text=element_text(size = 13, face = "plain"))
 ggsave(filename = "./Figure_Panels/EDFig.12d.pdf", width = 5.5, height = 5, useDingbats=FALSE)
 
