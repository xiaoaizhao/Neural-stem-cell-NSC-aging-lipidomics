## PCA plotting on all detected lipids in Primary NSC culture #2 dataset
rm(list=ls())
library(tidyverse)
library("scales")

# setwd(rstudioapi::getActiveProject())

########Import color LUT for all future plots with KO #######################################################################
load("./Output_Data/KO_LUT_paper.order2024-02-08.Rdata")

######PCA on individual lipids###########################################################################################
load("./Output_Data/Exp2_Norm_Impt_log2_all693_lipids.Rdata")

my.lipids <-  as.data.frame(t(Exp2_Impt_norm_conc_no_conc_all))
df_pca <- prcomp(my.lipids, center = TRUE, scale. = T)

df_out <- as.data.frame(df_pca$x)

df_out$group <- ""
df_out$group[grepl("N", rownames(df_out))] <- "Control"
df_out$group[grepl("A", rownames(df_out))] <- "Agpat3"
df_out$group[grepl("E", rownames(df_out))] <- "Elovl5"
df_out$group[grepl("F", rownames(df_out))] <- "Fads2"
df_out$group[grepl("M", rownames(df_out))] <- "Mboat2"
df_out$group[grepl("P", rownames(df_out))] <- "Pla2g4e"
df_out$age <- ""
df_out$age <- ifelse(grepl("Y", rownames(df_out)), "Young", "Old")
head(df_out)

xplot <- "PC1"
yplot <- "PC2"

df_out$age <- factor(df_out$age, levels = c("Young", "Old"))
df_out$group <- factor(df_out$group, levels = c("Control", "Elovl5", "Mboat2", "Fads2", "Agpat3", "Pla2g4e"))

p<-ggplot(df_out,aes(x=df_out[,xplot],y=df_out[,yplot],color= group))
p+geom_point(aes(shape = factor(age)), size=4, alpha=0.8)+
  theme_classic()+
  xlab(paste0(xplot, ": ", format(summary(df_pca)$importance[2,xplot] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0(yplot, ": ", format(summary(df_pca)$importance[2,yplot] * 100,
                                  digits = 3), " % variance"))+
  scale_color_manual(values = KO.LUT.paperorder$Color)+
  scale_shape_manual(values =c(16, 17)) +
  labs(title = "PCA Primary culture #2")+
  theme(text=element_text(size = 13, face = "plain"))
ggsave(filename = "./Figure_Panels/Fig.4c.pdf", width = 5, height = 5, useDingbats=FALSE)

####==== PCA on Mboat2 KO====
rm(list = ls())
pal4 <- c("cyan3", "magenta3")
load("./Output_Data/Exp2_Norm_Impt_log2_all693_lipids.Rdata")

pko <- c("_N", "_M")
mtx <- Exp2_Impt_norm_conc_no_conc_all %>% 
  dplyr::select(matches(pko))

my.lipids <-  as.data.frame(t(mtx))
df_pca <- prcomp(my.lipids, center = TRUE, scale. = T)

df_out <- as.data.frame(df_pca$x) %>% 
  mutate(., group = case_when(
    grepl("N", rownames(.)) ~ "Control",
    grepl("M", rownames(.)) ~ "Mboat2"
  )) %>% 
  mutate(., age = ifelse(grepl("Y", rownames(.)), "Young", "Old"))

xplot <- "PC1"
yplot <- "PC2"

df_out$age <- factor(df_out$age, levels = c("Young", "Old"))
df_out$group <- factor(df_out$group, levels = c("Control", "Mboat2"))

df_out <- df_out %>% 
  mutate(Cat = ifelse(group == "Control", "Control", "KO"))

p<-ggplot(df_out,aes(x=df_out[,xplot],y=df_out[,yplot]))
p+
  geom_point(size=3, alpha=0.8, stroke = 2, aes(color = age, shape = group))+
  theme_classic()+
  xlab(paste0(xplot, ": ", format(summary(df_pca)$importance[2,xplot] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0(yplot, ": ", format(summary(df_pca)$importance[2,yplot] * 100,
                                  digits = 3), " % variance"))+
  scale_color_manual(values = pal4) +
  scale_shape_manual(values = c(2, 11)) +
  labs(title = paste0(as.character(unique(df_out$group)), collapse = " "))+
  theme(text=element_text(size = 13, face = "plain")) +
  theme(axis.text=element_text(colour="black"))
ggsave(filename = "./Figure_Panels/Fig.4f.pdf", width = 5, height = 5, useDingbats=FALSE)