## Figure1  plot
### use z score for all lipids on Exp1, Exp2 and Batch 3 data
### Z scores were calculated by previous script - `PCA_Exp1_Exp2_Batch3_getZscore.R`

rm(list = ls())
load("./Output_Data/zscore.lpd.Exp2CTRL.for.cmb.PCA.Rdata")
load("./Output_Data/zscore.lpd.Exp3.for.cmb.PCA.Rdata")
load("./Output_Data/zscore.lpd.Exp1.for.cmb.PCA.Rdata")

####Subset to get common lipids across 3 experiments (96 lipids total) and generate a PCA#########################################################

PCA.3.datasets <- list( z.ko, z.lc, z.Exp3) %>% 
  reduce(inner_join, by = "LipidIon") %>%  #96 common lipids
  column_to_rownames(var = "LipidIon")

my.lipids <-  as.data.frame(t(PCA.3.datasets))
df_pca <- prcomp(my.lipids)

df_out <- as.data.frame(df_pca$x)

df_out <- df_out %>% 
  mutate(., 
         Cell_Type = ifelse(grepl("_N|qNSC|_Quiescent", row.names(df_out), ignore.case = T), "Quiescent", "Activated")) %>% 
  mutate(., Age = case_when(
    grepl("^Y", rownames(.)) ~ "Young",
    grepl("^O", rownames(.)) ~ "Old")) %>% 
  mutate(., Exp = case_when(
    grepl("_N", rownames(.)) ~ "Primary Culture #2",
    grepl("_Quiescent|_Activated", rownames(.)) ~ "Primary Culture #1",
    grepl("-A|-Q", rownames(.)) ~ "Primary Culture #3"
  )) %>% 
  relocate(c(Cell_Type, Age, Exp), .before = "PC1")

head(df_out)


xplot <- "PC1"
yplot <- "PC2"

df_out$Age <- factor(df_out$Age, levels = c("Young", "Old"))
df_out$Cell_Type <- factor(df_out$Cell_Type, levels = c("Activated", "Quiescent"))
df_out$Exp <- factor(df_out$Exp, levels = c("Primary Culture #1", "Primary Culture #2", "Primary Culture #3"))
df_out <- df_out %>% 
  mutate(., Exp_celltype = paste(Exp, Cell_Type, sep = "_"))

df_out$Exp_celltype <- factor(df_out$Exp_celltype, levels = c("Primary Culture #1_Activated",
                                                              "Primary Culture #1_Quiescent",
                                                              "Primary Culture #2_Quiescent",
                                                              "Primary Culture #3_Activated",
                                                              "Primary Culture #3_Quiescent"))

pal4 <- c("cyan3", "magenta3")
p<-ggplot(df_out,aes(x=df_out[,xplot],y=df_out[,yplot]))
p+geom_point(aes(shape = Exp_celltype, color= Age), size=3, alpha=0.8, fill = "grey70", stroke = 2)+
  theme_classic()+
  xlab(paste0(xplot, ": ", format(summary(df_pca)$importance[2,xplot] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0(yplot, ": ", format(summary(df_pca)$importance[2,yplot] * 100,
                                  digits = 3), " % variance"))+
  scale_color_manual(values = pal4)+
  scale_shape_manual(values =c(21,24,2,16,17))+
  theme(axis.text = element_text(colour = "black"))
# 
ggsave("./Figure_Panels/EDFig.1e.pdf", width = 7, height = 5,
       useDingbats=FALSE)

### ====== PCA plot on quiescent cells only ==================================================================
### ====== PCA plot on quiescent cells only ==================================================================
### ====== PCA plot on quiescent cells only ==================================================================
Q.3sets <- list(z.ko, z.Exp3, z.lc) %>% 
  reduce(inner_join, by = "LipidIon") %>%  
  column_to_rownames(var = "LipidIon") %>% 
  dplyr::select(-matches("Activated|aNSC"))

my.lipids <-  as.data.frame(t(Q.3sets))
df_pca <- prcomp(my.lipids)

df_out <- as.data.frame(df_pca$x)

df_out<- df_out %>% 
  mutate(., Age = case_when(
    grepl("^Y", rownames(.)) ~ "Young",
    grepl("^O", rownames(.)) ~ "Old")) %>% 
  mutate(., Exp = case_when(
    grepl("_Quiescent", rownames(.)) ~ "Primary Culture #1",
    grepl("_N", rownames(.)) ~ "Primary Culture #2",
    grepl("_qNSC-Q", rownames(.)) ~ "Primary Culture #3"
  )) %>% 
  relocate(c(Age, Exp), .before = "PC1")

df_out$Age <- factor(df_out$Age, levels = c("Young", "Old"))

head(df_out)

xplot <- "PC1"
yplot <- "PC2"

df_out$Exp <- factor(df_out$Exp, levels = c(
  "Primary Culture #1", 
  "Primary Culture #2", 
  "Primary Culture #3"))


p<-ggplot(df_out,aes(x=df_out[,xplot],y=df_out[,yplot]))
p+geom_point(aes(color= Age, shape = Exp), size=3, alpha=0.8, fill = "grey70", stroke = 2)+
  theme_classic()+
  xlab(paste0(xplot, ": ", format(summary(df_pca)$importance[2,xplot] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0(yplot, ": ", format(summary(df_pca)$importance[2,yplot] * 100,
                                  digits = 3), " % variance"))+
  scale_color_manual(values = pal4)+
  scale_shape_manual(values =c(24, 2, 17)) 

ggsave("./Figure_Panels/EDFig.1f.pdf", width = 6, height = 5,
       useDingbats=FALSE)
