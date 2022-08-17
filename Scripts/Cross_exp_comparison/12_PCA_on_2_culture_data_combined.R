##Project Primary Culture #1 dataset and Primary Culture #2 dataset into one PCA space

#Steps:
#1. Import log2 transformed data from both experiment
#2. Calculate Z score on all lipids for each experiment
#3. Subset to get common lipids across both experiment (266 lipids total) and generate a PCA

rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

####Load log2 transformed data from 1st Primary NSC culture data####################################################################################
load(file = "./Output_data/Spike-in_norm_MedNorm_all_373_lipid.Rdata") #this data is already log2 transformed

##calculate z score for Primary Culture #1 dataset
z.lc <- Impt_norm_373all %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Samples", values_to = "Int") %>%
  group_by(LipidIon) %>% 
  mutate(z = as.vector(scale(Int))) %>% 
  pivot_wider(-Int, names_from = Samples, values_from = z) 

####Load log2 transformed data from Primary Culture #2 dataset####################################################################################
load(file = "./Output_Data/Exp2_Norm_Impt_log2_all693_lipids.Rdata")
df.ctrl <- Exp2_Impt_norm_conc_no_conc_all %>%
  select(., contains("_N"))

##calculate z score for Primary Culture #2 dataset, control samples only
z.ko <- df.ctrl %>%
  rownames_to_column(., var = "LipidIon") %>%
  pivot_longer(-LipidIon, names_to = "Samples", values_to = "Int") %>%
  group_by(LipidIon) %>% 
  mutate(z = as.vector(scale(Int))) %>% 
  pivot_wider(-Int, names_from = Samples, values_from = z) 


####Subset to get common lipids across both experiment (266 lipids total) and generate a PCA#########################################################
z.ovlp <- inner_join(z.ko, z.lc, by = "LipidIon") %>%
  column_to_rownames(., var = "LipidIon") ##266 common lipids between 2 datasets

my.lipids <-  as.data.frame(t(z.ovlp))
df_pca <- prcomp(my.lipids)

df_out <- as.data.frame(df_pca$x)

df_out$`Cell Type` <- ""
df_out$`Cell Type` <- ifelse(grepl("qui|_N", row.names(df_out), ignore.case = T), "Quiescent", "Activated")
df_out$Age <- ""
df_out$Age <- ifelse(grepl("Y", rownames(df_out)), "Young", "Old")
df_out$Exp <- ifelse(grepl("_N", rownames(df_out)), "KO", "LC-MS")
head(df_out)

xplot <- "PC1"
yplot <- "PC2"

df_out$Age <- factor(df_out$Age, levels = c("Young", "Old"))
df_out$`Cell Type` <- factor(df_out$`Cell Type`, levels = c("Activated", "Quiescent"))
df_out$Exp <- factor(df_out$Exp, levels = c("LC-MS", "KO"))
df_out <- df_out %>%
  mutate(., Exp_celltype = paste(Exp, `Cell Type`, sep = "_")) %>%
  mutate(., Exp_celltype = case_when(
    grepl("KO", Exp_celltype) ~ "Primary Culture #2 - qNSC", 
    grepl("LC-MS_Activated", Exp_celltype) ~ "Primary Culture - aNSC",
    grepl("LC-MS_Quiescent", Exp_celltype) ~ "Primary Culture - qNSC",
  ))

df_out$Exp_celltype <- factor(df_out$Exp_celltype, levels = c("Primary Culture #2 - qNSC",
                                                              "Primary Culture - aNSC",
                                                              "Primary Culture - qNSC"))

palette1 <- c("darkgoldenrod", "maroon")

p<-ggplot(df_out,aes(x=df_out[,xplot],y=df_out[,yplot]))
p+geom_point(aes(shape = `Cell Type`, color= Age), size=4, alpha=0.8)+
  theme_classic()+
  xlab(paste0(xplot, ": ", format(summary(df_pca)$importance[2,xplot] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0(yplot, ": ", format(summary(df_pca)$importance[2,yplot] * 100,
                                  digits = 3), " % variance"))+
  scale_colour_manual(values = palette1)+
  scale_shape_manual(values =c(16, 17)) +
  stat_ellipse(aes(linetype = Exp_celltype))
ggsave(paste0("./Figure_Panels/Fig_S1c.pdf"), width = 6, height = 5,
       useDingbats=FALSE)
