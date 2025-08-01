## PCA analysis on lipidomics of qNSC Mboat2 overexpression in vitro

rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")

load("./Output_data/M2PM_Norm_Impt_log2_conc_275_noconc_48_lipids.Rdata")

Smp.Key <- read.csv("~/Dropbox/Stanford/Lipidomics/2023_Feb_lipidomics/March_Sample_list_071123_forR.csv", stringsAsFactors = F)
Smp.Key.e <- Smp.Key %>% 
  select(Sample.Name, Sample_ID)  %>% 
  rowwise() %>% 
  mutate(Exp = case_when(
    grepl("_EGFP|_Mb2_OE", Sample.Name) ~ "Mboat2_OE",
    grepl("_meth|lpd", Sample.Name) ~ "GPMV_Sup",
  ))  

M2OE.lg2 <- M2PM.Impt_norm_conc_no_conc_all %>% 
  rename_at(vars(matches(Smp.Key.e$Sample_ID)), ~Smp.Key.e$Sample.Name) %>% 
  select(matches("_EGFP|_Mb2_OE"))

my.lipids <-  as.data.frame(t(M2OE.lg2))
df_pca <- prcomp(my.lipids, center = TRUE, scale. = T)

df_out <- as.data.frame(df_pca$x) %>% 
  rownames_to_column(., var = "Lpd_smple")

df_out.p <- df_out %>% 
  rowwise() %>% 
  mutate(OE_cond = case_when(
    grepl("_Mb2_OE", Lpd_smple) ~ "Mboat2 OE",
    grepl("_EGFP", Lpd_smple) ~ "Control"
  ))  %>%
  mutate(., Culture = substr(Lpd_smple, 1, str_locate(Lpd_smple, "_")-1)) %>% 
  mutate(., Age = ifelse(grepl("^Y", Lpd_smple), "Young", "Old")) %>% 
  relocate(c(Culture, Age, OE_cond, Lpd_smple), .after = Lpd_smple)

xplot <- "PC1"
yplot <- "PC3"

df_out.p <- df_out.p %>% 
  ungroup() 

df_out.p$Age <- factor(df_out.p$Age, levels = c("Young", "Old"))
df_out.p$OE_cond <- factor(df_out.p$OE_cond, levels = c("Control", "Mboat2 OE"))
df_out.p <- as.data.frame(df_out.p)

pal4 <- c("cyan3", "magenta3")

p<-ggplot(df_out.p,aes(x=.data[[xplot]],y=.data[[yplot]]))
p+geom_point(size=3, alpha=0.8, stroke = 2, aes(color = Age, shape = OE_cond))+
  theme_classic()+
  xlab(paste0(xplot, ": ", format(summary(df_pca)$importance[2,xplot] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0(yplot, ": ", format(summary(df_pca)$importance[2,yplot] * 100,
                                  digits = 3), " % variance"))+
  scale_color_manual(values = pal4) +
  scale_shape_manual(values = c(5, 7)) +
  labs(title = "Mboat2 OE")+
  theme(text=element_text(size = 13, face = "plain"))
ggsave(filename =paste0("./Figure_Panels/Fig.4j.pdf"), width = 5.2, height = 5, useDingbats=FALSE)