
##PCA following spike-in normalization, median concentration normalization and imputation
rm(list=ls())
library(tidyverse)
library(pROC) 
library(data.table)
setwd(rstudioapi::getActiveProject())

######PCA #############################################################################################
##130 lipids total 

load(file = "./Output_Data/Invivo_Norm_Impt_log2_all130_lipids.Rdata") #this data is already log2 transformed
my.lipids <-  as.data.frame(t(Impt_norm_conc_no_conc_all))
df_pca <- prcomp(my.lipids, center = TRUE, scale. = T)

df_out <- as.data.frame(df_pca$x)

df_out$age <- ""
df_out$age <- ifelse(grepl("Y", rownames(df_out)), "Young", "Old")
df_out$age <- factor(df_out$age, levels = c("Young", "Old"))
head(df_out)

mDF <- df_out %>% 
  select(., -age) %>% 
  rownames_to_column(., var = "Samples") %>% 
  pivot_longer(-Samples, names_to = "PC", values_to = "Loading") %>% 
  mutate(., age = ifelse(grepl("O", Samples), "Old", "Young"))

mDF$age <- factor(mDF$age, levels = c("Young", "Old"))

mDFt <- mDF %>% 
  group_by(PC) %>%
  summarize(Wilcox_Pval = wilcox.test(Loading ~ age)$p.value) %>%
  mutate(., Padj = p.adjust(Wilcox_Pval, "BH")) %>% 
  arrange(., Wilcox_Pval)

AUC.all <- list()

for (PCn in mDF$PC) {
  AUC.all[[PCn]] <- mDF %>% 
    filter(., PC == PCn) %>%
    summarise(., AUC = as.numeric(roc(age, Loading)$auc)) %>% 
    mutate(., PCname = PCn)
}

AUC <- bind_rows(AUC.all) %>% 
  arrange(desc(AUC))
  
AUC

xplot = "PC6"
yplot = "PC9"

df.tb <- as.data.table(df_out)
decB <- glm(age ~ df.tb[, get(xplot)] + df.tb[, get(yplot)],
            data   = df.tb,
            family = binomial)

summary(decB)

slope <- coef(decB)[2]/(-coef(decB)[3])
intercept <- coef(decB)[1]/(-coef(decB)[3]) 


a<- levels(factor(df_out$group))

cbp2 <- c("darkgoldenrod", "maroon")

p<-ggplot(df_out,aes(x=df_out[,xplot],y=df_out[,yplot],color= age))
p+geom_point( size=4, alpha=0.8)+
  theme_classic()+
  xlab(paste0(xplot, ": ", format(summary(df_pca)$importance[2,xplot] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0(yplot, ": ", format(summary(df_pca)$importance[2,yplot] * 100,
                                  digits = 3), " % variance"))+
  scale_colour_manual(values = cbp2)+
  scale_shape_manual(values =c(16, 17)) +
  labs(title = "In vivo sorted qNSC", color = "Age")+
  geom_abline(slope = slope, intercept = intercept, linetype = 2)+
  theme(text=element_text(size = 13, face = "plain"),
        axis.text = element_text(colour = "black"))
ggsave(paste0("./Figure_Panels/Fig_S2b.pdf"), width = 5, height = 5, useDingbats=FALSE)


