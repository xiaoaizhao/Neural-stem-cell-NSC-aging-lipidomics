### PCA between Primary culture #3 nad Lipidyzer qNSCs

setwd(rstudioapi::getActiveProject())
rm(list=ls())
library(tidyverse)

## Quiescent cells only
### ====quiescent cells only - lipidyzer ====####
load("./Output_Data/Lipidyzer_qNSC_for.LC-MS.ovlp.Rdata") #this data is NOT log2 transformed
#### log 2 transformation
D.qui.log2 <- Ldz.Qui.nodup %>% 
  mutate(across(matches("^Y|^O", ignore.case = FALSE), log2))

####calculate z score for quiescent cells of Lipidyzer

z.q.dyzer <- D.qui.log2 %>%
  pivot_longer(-c(LipidIon, ID_string), names_to = "Samples", values_to = "Conc") %>%
  group_by(LipidIon) %>% 
  mutate(z = as.vector(scale(Conc))) %>% 
  pivot_wider(-Conc, names_from = Samples, values_from = z) 

####====quiescent cells only - Batch 3 ====####
load("./Output_Data/Exp3_lipid.conc.reformat.for.lipidyzer.overlap.Rdata")

#### log 2 transformation
E3.qui.log2 <- E3.lpd.conc.fmt %>% 
  select(-LipidIon1) %>% 
  select(matches("LipidIon|ID_string|_qNSC-Q")) %>% 
  mutate(across(matches("^Y|^O", ignore.case = FALSE), log2))

####calculate z score for quiescent cells of Batch 3

z.q.E3 <- E3.qui.log2 %>%
  pivot_longer(-c(LipidIon, ID_string), names_to = "Samples", values_to = "Conc") %>%
  group_by(LipidIon) %>% 
  mutate(z = as.vector(scale(Conc))) %>% 
  pivot_wider(-Conc, names_from = Samples, values_from = z) 

### ====== PCA plot on quiescent cells only ==================================================================
Q.2sets <- list(z.q.E3, z.q.dyzer) %>% 
  reduce(inner_join, by = "ID_string") %>%  #103 common lipids
  column_to_rownames(var = "LipidIon.x") %>% 
  dplyr::select(-matches("LipidIon|ID_string"))

my.lipids <-  as.data.frame(t(Q.2sets))
df_pca <- prcomp(my.lipids)

df_out <- as.data.frame(df_pca$x)

df_out<- df_out %>% 
  mutate(., Age = case_when(
    grepl("^Y", rownames(.)) ~ "Young",
    grepl("^O", rownames(.)) ~ "Old")) %>% 
  mutate(., Exp = case_when(
    grepl("_qui", rownames(.)) ~ "Lipidyzer",
    grepl("_qNSC-Q", rownames(.)) ~ "Primary Culture #3"
  )) %>% 
  relocate(c(Age, Exp), .before = "PC1")

df_out$Age <- factor(df_out$Age, levels = c("Young", "Old"))

head(df_out)

xplot <- "PC1"
yplot <- "PC2"

df_out$Exp <- factor(df_out$Exp, levels = c(
  "Primary Culture #3", 
  "Lipidyzer"))

pal4 <- c("cyan3", "magenta3")

p<-ggplot(df_out,aes(x=df_out[,xplot],y=df_out[,yplot]))
p+geom_point(aes(color= Age, shape = Exp), size=3, alpha=0.8, fill = "grey70", stroke = 2)+
  theme_classic()+
  xlab(paste0(xplot, ": ", format(summary(df_pca)$importance[2,xplot] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0(yplot, ": ", format(summary(df_pca)$importance[2,yplot] * 100,
                                  digits = 3), " % variance"))+
  scale_color_manual(values = pal4)+
  scale_shape_manual(values =c(17, 10)) +
  theme(axis.text = element_text(colour = "black"))

ggsave(filename = "./Figure_Panels/EDFig.3f.pdf", width = 6, height = 5,
       useDingbats=FALSE)


