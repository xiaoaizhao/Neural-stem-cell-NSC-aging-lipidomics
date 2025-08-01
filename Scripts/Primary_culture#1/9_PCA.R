##PCA following spike-in normalization and median intensity normalization
##Result: PC1 separates cell type, PC3 separates age.
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

######PCA #############################################################################################
##373
load(file = "./Output_data/Spike-in_norm_MedNorm_all_373_lipid.Rdata") #this data is already log2 transformed

my.lipids <-  as.data.frame(t(Impt_norm_373all))
df_pca <- prcomp(my.lipids, center = TRUE, scale. = T)

df_out <- as.data.frame(df_pca$x)

df_out$group <- ""
df_out$group <- ifelse(grepl("qui", row.names(df_out), ignore.case = T),"Quiescent", "Activated")
df_out$group <- factor(df_out$group, levels = c("Activated", "Quiescent"))
df_out$age <- ""
df_out$age <- ifelse(grepl("Y", rownames(df_out)), "Young", "Old")
df_out$ExpDate <- ifelse(grepl("08", rownames(df_out)), "08", "11")
df_out$age <- factor(df_out$age, levels = c("Young", "Old"))
df_out$group <- factor(df_out$group, levels = c("Activated", "Quiescent"))
head(df_out)

####Plotting#####
xplot <- "PC1"
yplot <- "PC3"
pal4 <- c("cyan3", "magenta3")
p<-ggplot(df_out,aes(x=df_out[,xplot],y=df_out[,yplot]))
p+geom_point(aes(shape = group, color= age), size=3, alpha=0.8, fill = "grey70", stroke = 1)+
  # geom_point(aes(shape = group)), size=4, alpha=0.8)+
  theme_classic()+
  xlab(paste0(xplot, ": ", format(summary(df_pca)$importance[2,xplot] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0(yplot, ": ", format(summary(df_pca)$importance[2,yplot] * 100,
                                  digits = 3), " % variance"))+
  scale_colour_manual(values = pal4)+
  scale_shape_manual(values =c(21, 24)) +
  labs(title = "Primary culture #1", color = "Age", shape = "Cell Type")+
  theme(text=element_text(size = 13, face = "plain"))
ggsave(filename = "./Figure_Panels/EDFig.2.PCA.Primary_culture_#1.pdf", width = 5, height = 5, useDingbats=FALSE)
