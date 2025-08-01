#' ### PCA on Experiment 2 data

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
#' 
#' Load in log2 transformed matrix of all samples
#' 
## ---------------------------------------------------------------------------------------------------------------------
load("./Output_Data/Exp2_Norm_Impt_log2_all693_lipids.Rdata")

#' 
#' #### Control cell on its own
## ---------------------------------------------------------------------------------------------------------------------
Ctrl <- Exp2_Impt_norm_conc_no_conc_all %>% 
  select(contains("_N"))

my.lipids <-  as.data.frame(t(Ctrl))
df_pca <- prcomp(my.lipids, center = TRUE, scale. = T)

df_out <- as.data.frame(df_pca$x)

df_out$age <- ""
df_out$age <- ifelse(grepl("Y", rownames(df_out)), "Young", "Old")
# head(df_out)

xplot <- "PC1"
yplot <- "PC2"

df_out$age <- factor(df_out$age, levels = c("Young", "Old"))

pal4 <- c("cyan3", "magenta3")

p<-ggplot(df_out,aes(x=df_out[,xplot],y=df_out[,yplot],color= age))
p+geom_point(size=3, alpha=0.8, shape = 2, stroke = 2)+
  theme_classic()+
  xlab(paste0(xplot, ": ", format(summary(df_pca)$importance[2,xplot] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0(yplot, ": ", format(summary(df_pca)$importance[2,yplot] * 100,
                                  digits = 3), " % variance"))+
  scale_color_manual(values = pal4)+
  labs(title = "Primary culture #2 - controls")+
  theme(text=element_text(size = 13, face = "plain"))
  ggsave(filename = "./Figure_Panels/EDFig.2.PCA.Primary_culture_#2_control.pdf", width = 5, height = 5, useDingbats=FALSE)

#' 
#' #### Control cell paired with each KO
## ---------------------------------------------------------------------------------------------------------------------
KO.ls <- c("_A", "_E", "_F", "_P", "_M")

for (ko.n in KO.ls) {
  pko <- c("_N", ko.n)
  mtx <- Exp2_Impt_norm_conc_no_conc_all %>% 
    dplyr::select(matches(pko))
  
  my.lipids <-  as.data.frame(t(mtx))
  df_pca <- prcomp(my.lipids, center = TRUE, scale. = T)

  df_out <- as.data.frame(df_pca$x) %>% 
    mutate(., group = case_when(
      grepl("N", rownames(.)) ~ "Control",
      grepl("A", rownames(.)) ~ "Agpat3",
      grepl("E", rownames(.)) ~ "Elovl5",
      grepl("F", rownames(.)) ~ "Fads2",
      grepl("M", rownames(.)) ~ "Mboat2",
      grepl("P", rownames(.)) ~ "Pla2g4e",
    )) %>% 
    mutate(., age = ifelse(grepl("Y", rownames(.)), "Young", "Old"))

    xplot <- "PC1"
    yplot <- "PC2"
    
    df_out$age <- factor(df_out$age, levels = c("Young", "Old"))
    df_out$group <- factor(df_out$group, levels = c("Control", "Elovl5", "Mboat2", "Fads2", "Agpat3", "Pla2g4e"))
    
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
    ggsave(filename = paste0("./Figure_Panels/EDFig.2.PCA.Primary_culture_#2_",paste0(as.character(unique(df_out$group)), collapse = "_"), ".pdf"), width = 5, height = 5, useDingbats=FALSE)
}

#' 
