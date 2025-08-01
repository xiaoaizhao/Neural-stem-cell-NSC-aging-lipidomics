## PCA analysis on Primary culture #3 samples
## Steps:
#1. Even with adjusted concentration of O8 sample (Script #6), it still appears as an outlier on the PCA. Will remove it for all downstream analysis
#2. PCA does not look amazing. Let's see how the effect size comparison will shake out.
library(tidyverse)
# library(ggthemes)
# library("scales")
rm(list=ls())

# Use log2 transformed data for PCA analysis
load("./Output_data/Exp3_Norm_Impt_log2_conc_356_noconc_60_lipids.Rdata")
  
Smp.Key <- read.csv("./Input_Data/Batch3_sample_key_forR.csv", stringsAsFactors = F)
Smp.Key.e <- Smp.Key %>% 
  select(Sample.Name, Sample_ID) %>% 
  filter(!grepl("PBS", Sample_ID))

## Annotate sample ID with sample name
Exp3.df <- Exp3.Impt_norm_conc_no_conc_all %>% 
  rename_at(vars(matches(Smp.Key.e$Sample_ID)), ~Smp.Key.e$Sample.Name) %>% 
  filter(!grepl("Cholesterol\\+H\\-H2O_23_Positive|Ch\\+H\\-H2O", rownames(.))) 

## sanity check
l=29
n=15
Exp3.df[l,n] == Exp3.Impt_norm_conc_no_conc_all[rownames(Exp3.Impt_norm_conc_no_conc_all) == rownames(Exp3.df)[l],
                                         Smp.Key.e$Sample_ID[Smp.Key.e$Sample.Name == colnames(Exp3.df)[n]]]

### O8-aNSC sample was accidentally diluted 4 times during sample preparation. Although concentration was adjusted, this sample still looks to be an outlier from all other samples and therefore is excluded from following analyses
## Remove O8 aNSC sample from all future analysis
Exp3.no.O8A <- Exp3.df %>% 
  select(-`O8_aNSC-A`)

save(Exp3.no.O8A, file = "./Output_Data/Exp3_raw_rmv_O3ANSC.Rdata")


## PCA analysis
## ====Quiescent and activated samples together====
Exp3.Q <- Exp3.no.O8A %>% 
  select(contains("qNSC-Q"))

my.lipids <-  as.data.frame(t(Exp3.no.O8A))

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
  ggtitle("Primary culture #3")
ggsave(filename = "./Figure_Panels/Fig.1b.pdf", width = 5.5, height = 5, useDingbats=FALSE)


## ====Quiescent samples only====
my.lipids <-  as.data.frame(t(Exp3.Q))

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

xplot <- "PC2"
yplot <- "PC3"

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
  ggtitle("Primary culture #3 qNSCs only") +
  theme(text=element_text(size = 13, face = "plain")) +
  scale_shape_manual(values =c(17)) +
  theme(axis.text=element_text(colour="black"))
  # stat_ellipse(aes(color = Age), type = "t", linetype = 2)
ggsave(filename = "./Figure_Panels/Fig.1c.pdf", width = 5.5, height = 5, useDingbats=FALSE)

