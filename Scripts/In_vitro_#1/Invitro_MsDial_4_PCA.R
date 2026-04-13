### PCA

library(tidyverse)
rm(list=ls())

# Use log2 transformed data for PCA analysis
load("./Output_data/Invitro_MsD.Norm_Impt_log2_conc_404_lipids.Rdata")

### O8-aNSC sample was accidentally diluted 4 times during sample preparation. Although concentration was adjusted, this sample still looks to be an outlier from all other samples and therefore is excluded from following analyses
## Remove O8 aNSC sample from all future analysis
Invitro.no.O8A.MsD <- Invitro.log2.impt.norm.conc.MsD %>%
  select(-`O8_aNSC-A`)
 
save(Invitro.no.O8A.MsD, file = "./Output_Data/Invitro.MsD.404.log2.norm.impt.no.O8aNSC.Rdata")


## PCA analysis
## ====Quiescent and activated samples together====
Invitro.Q <- Invitro.no.O8A.MsD %>% 
  select(contains("qNSC-Q"))

my.lipids <-  as.data.frame(t(Invitro.no.O8A.MsD))

df_pca <- prcomp(my.lipids, center = TRUE, scale. = T)

df_out <- as.data.frame(df_pca$x) %>% 
  rownames_to_column(., var = "Key")

df_out.p <- df_out %>% 
  rowwise() %>% 
  mutate(., Age = ifelse(grepl("^Y", Key), "Young", "Old")) %>% 
  mutate(., CellType = case_when(
    grepl("qNSC-Q", Key) ~ "qNSC",
    grepl("aNSC-A", Key) ~ "aNSC"
  )) %>% 
  relocate(c(Age, CellType), .after = Key)

xplot <- "PC1"
yplot <- "PC2"

df_out.p$Age <- factor(df_out.p$Age, levels = c("Young", "Old"))
df_out.p$CellType <- factor(df_out.p$CellType, levels = c("aNSC", "qNSC"))

pal4 <- c("cyan3", "magenta3")
p<-ggplot(df_out.p,aes(x=.data[[xplot]],y=.data[[yplot]]
                       # ,label = Key
                       ))
p+geom_point(
  aes(shape = CellType, color = Age),
  size=4, alpha=0.8)+
  theme_classic()+
  xlab(paste0(xplot, ": ", format(summary(df_pca)$importance[2,xplot] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0(yplot, ": ", format(summary(df_pca)$importance[2,yplot] * 100,
                                  digits = 3), " % variance"))+
  scale_shape_manual(values =c(16, 17)) +
  scale_color_manual(values = pal4)+
  labs(shape = "Cell Type", color = "Age")+
  theme(text=element_text(size = 13, face = "plain"))+
  ggtitle("In vitro") 
ggsave(filename = "./Figure_Panels/Fig.1b.pdf", width = 5.5, height = 5, useDingbats=FALSE)


## ====Quiescent samples only====
my.lipids <-  as.data.frame(t(Invitro.Q))

df_pca <- prcomp(my.lipids, center = TRUE, scale. = T)

df_out <- as.data.frame(df_pca$x) %>% 
  rownames_to_column(., var = "Key")

df_out.p <- df_out %>% 
  rowwise() %>% 
  mutate(., Age = ifelse(grepl("^Y", Key), "Young", "Old")) %>% 
  mutate(., CellType = case_when(
    grepl("qNSC-Q", Key) ~ "qNSC",
    grepl("aNSC-A", Key) ~ "aNSC"
  )) %>% 
  relocate(c(Age, CellType), .after = Key)

xplot <- "PC1"
yplot <- "PC2"

df_out.p$Age <- factor(df_out.p$Age, levels = c("Young", "Old"))
df_out.p$CellType <- factor(df_out.p$CellType, levels = c("aNSC", "qNSC"))


p<-ggplot(df_out.p,aes(x=.data[[xplot]],y=.data[[yplot]]
                       # ,label = Key
))
p+geom_point(
  aes(color = Age, shape = CellType),
  size=4, alpha=0.8)+
  theme_classic()+
  xlab(paste0(xplot, ": ", format(summary(df_pca)$importance[2,xplot] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0(yplot, ": ", format(summary(df_pca)$importance[2,yplot] * 100,
                                  digits = 3), " % variance"))+
  scale_color_manual(values = pal4)+
  ggtitle("In vitro qNSCs only") +
  theme(text=element_text(size = 13, face = "plain")) +
  scale_shape_manual(values =c(17)) +
  theme(axis.text=element_text(colour="black"))
  
ggsave(filename = "./Figure_Panels/Fig.1c.pdf", width = 5.5, height = 5, useDingbats=FALSE)

