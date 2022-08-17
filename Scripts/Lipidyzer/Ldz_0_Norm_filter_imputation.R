##Lipidyzer data on primary NSC culture
rm(list=ls())
library(tidyverse)

setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Pre-processing_functions.R")

#### Append sample name to data frame ####===============================================================
# Upload data
all <- read.csv("./Input_Data/a_L_species_Conc_BlkNorm.csv",
                stringsAsFactors = F) 
#exclude QC + Blk, import DNA normalization factor
all.df <- all %>%
  select(., -matches("QC|B"))
#import sample name spreadsheet
key <- read.csv("./Input_Data/sample_key.csv", stringsAsFactors = F)
key <- key %>%
  mutate(., ID = paste0("X", Key))
#Append sample name to data matrix
all_features <- all.df %>%
  rename_at(vars(key$ID), ~key$Samples) %>%
  column_to_rownames(., var = "features")

save(all_features, file = paste0("./Output_Data/Lipidyzer_all_feature.Rdata"))

#### Median concentration normalization ####=============================================================
load("./Output_Data/Lipidyzer_all_feature.Rdata")
Med_ind <- apply(all_features, 2, median, na.rm=T)
Med_norm_factor <- as.numeric(scale(Med_ind, center = F, scale = median(Med_ind)))

Med_norm.ldz <- sweep(all_features, 2, Med_norm_factor, "/")
N_med <- apply(Med_norm.ldz, 2, median, na.rm=T)

save(Med_norm.ldz, file = paste0("./Output_Data/Ldz_Med_conc_normed_all_features.Rdata"))

#### Filter and Imputation####===========================================================================
##filter lipid should be detected at least in half of the samples (NA <10)
load("./Output_Data/Ldz_Med_conc_normed_all_features.Rdata")
Norm_filter <- filter.Ldz(Med_norm.ldz, 21) ##511 lipids meet this criteria

##Imputation to replace missing value
Norm_filter.ipt <- impt.Ldz(Norm_filter)

##Check the number of values that were replace by imputation = number of missing value to begin with
impt.check.Ldz(Norm_filter, Norm_filter.ipt)

Ldz.norm.impt <- Norm_filter.ipt
save(Ldz.norm.impt, file = paste0("./Output_Data/Ldz_511_lipids_log2_MedConc_norm_imputed.Rdata"))

##Also create a dataframe with values transform back to raw concentration####
ldz.raw_conc <- 2^(Ldz.norm.impt)
save(ldz.raw_conc, file = paste0("./Output_data/Ldz_511_lipids_backtoraw.conc_MedConc_norm_imputed.Rdata"))


####Plot PCA####===========================================================================================
load("./Output_Data/Ldz_511_lipids_log2_MedConc_norm_imputed.Rdata")
my.lipids <-  as.data.frame(t(Ldz.norm.impt))
df_pca <- prcomp(my.lipids, center = TRUE, scale. = T)

df_out <- as.data.frame(df_pca$x)

df_out$group <- ""
df_out$group <- ifelse(grepl("qui", row.names(df_out), ignore.case = T), "Quiescent", "Activated")
df_out$group <- factor(df_out$group, levels = c("Activated", "Quiescent"))
df_out$age <- ""
df_out$age <- ifelse(grepl("Y", rownames(df_out)), "Young", "Old")
head(df_out)

xplot <- "PC1"
yplot <- "PC2"
Palette1 <- c("darkgoldenrod", "maroon")
df_out$age <- factor(df_out$age, levels = c("Young", "Old"))

p<-ggplot(df_out,aes(x=df_out[,xplot],y=df_out[,yplot],color= age))
p+geom_point(aes(shape = factor(group)), size=4, alpha=0.8)+
  theme_classic()+
  xlab(paste0(xplot, ": ", format(summary(df_pca)$importance[2,xplot] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0(yplot, ": ", format(summary(df_pca)$importance[2,yplot] * 100,
                                  digits = 3), " % variance"))+
  scale_colour_manual(values = Palette1)+
  scale_shape_manual(values =c(16, 17)) +
  labs(title = "", color = "Age", shape = "Cell Type")+
  theme(text=element_text(size = 12))
ggsave(paste0("./Figure_Panels/Fig_S1f.pdf"), width = 5, height = 5,
       useDingbats=FALSE)

