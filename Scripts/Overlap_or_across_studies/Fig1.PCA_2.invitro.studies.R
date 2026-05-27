### PCA plot on In vitro and In vitro Experiment #2 together

rm(list = ls())
library(tidyverse)
source("./Scripts/Function_scripts/Pre-processing_functions.R")
load("./Output_Data/zscore.lpd.Exp2CTRL.for.cmb.PCA.MsD.frmt.Rdata")
load("./Output_Data/zscore.lpd.Invitro.for.cmb.PCAMsD.frmt.Rdata")
z.Invitro <- z.Invitro %>% 
  ungroup()
z.Invitro.rmvabc <- MsD.lpd.format(z.Invitro)
z.ko.ctrl <- z.ko %>% 
  ungroup()
z.exp2.rmvabc <- MsD.lpd.format(z.ko.ctrl)
####PCA on common lipids between Exp2 and quiescent cells of Exp3##########################################

PCA.2.datasets <- list( z.Invitro.rmvabc, z.exp2.rmvabc) %>% 
  reduce(inner_join, by = "LipidIon") %>%  #142 common lipids
  column_to_rownames(var = "LipidIon")

my.lipids <-  as.data.frame(t(PCA.2.datasets))
df_pca <- prcomp(my.lipids)

df_out <- as.data.frame(df_pca$x)

df_out<- df_out %>% 
  mutate(., Age = case_when(
    grepl("^Y", rownames(.)) ~ "Young",
    grepl("^O", rownames(.)) ~ "Old")) %>% 
  mutate(., Exp = case_when(
    grepl("_N", rownames(.)) ~ "In vitro Experiment #2",
    grepl("_qNSC-Q", rownames(.)) ~ "In vitro"
  )) %>% 
  relocate(c(Age, Exp), .before = "PC1")

df_out$Age <- factor(df_out$Age, levels = c("Young", "Old"))

head(df_out)

xplot <- "PC1"
yplot <- "PC2"
pal4 <- c("cyan3", "magenta3")


df_out$Exp <- factor(df_out$Exp, levels = c(
  "In vitro", 
  "In vitro Experiment #2"))


p<-ggplot(df_out,aes(x=df_out[,xplot],y=df_out[,yplot]))
p+geom_point(aes(color= Age, shape = Exp), size=3, alpha=0.8, fill = "grey70", stroke = 2)+
  theme_classic()+
  xlab(paste0(xplot, ": ", format(summary(df_pca)$importance[2,xplot] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0(yplot, ": ", format(summary(df_pca)$importance[2,yplot] * 100,
                                  digits = 3), " % variance"))+
  scale_color_manual(values = pal4)+
  scale_shape_manual(values =c(17,2)) +
  labs(title = "2 In vitro studies")

ggsave("./Figure_Panels/fig.S1E.pdf", width = 5.5, height = 5,
       useDingbats=FALSE)
