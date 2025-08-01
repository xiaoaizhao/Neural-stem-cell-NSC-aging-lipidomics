## Quantification of phospholipid order index from GPMV lipidomics
## Quantification of cholesterol concentration from GPMV lipidomics

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)

####====Phospholipid order index====
load("./Output_Data/GPMV_ClassSum_565_lipids.Rdata")
sum.df <- GPMV.classsum %>%
  pivot_wider(., names_from = "Class", values_from = "Class_sum") %>%
  mutate(., POI = (PE+SM)/PC) %>%
  mutate(., Age = ifelse(grepl("Y", Sample), "Young", "Old"))

sum.df$Age <- factor(sum.df$Age, levels = c("Young", "Old"))

pal4 <- c("cyan3", "magenta3")
c <- ggplot(sum.df, aes(Age, POI, color = Age))
c+  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.1, size = 3, alpha = 0.7)+
  theme_classic()+
  stat_compare_means(aes(group = Age), label = "p.format")+
  labs(title = "GPMV Phospho Order Index", y = "(PE+SM)/PC")+
  scale_color_manual(values = pal4)+
  theme(legend.position= "none")+
  theme(text=element_text(size = 13, face = "plain"))
ggsave(paste0("./Figure_Panels/Fig.3e.pdf"), width = 3.5, height = 5, useDingbats=FALSE)

sum.df.med <- sum.df %>% 
  group_by(Age) %>% 
  summarise(., MedianPOI = median(POI))

sum.df.med

####====Cholesterol concentration====
chol.GPMV <- sum.df %>% 
  select(Sample, Cholesterol, Age)

chol.GPMV$Age <- factor(chol.GPMV$Age, levels = c("Young", "Old"))

pal4 <- c("cyan3", "magenta3")

d <- ggplot(chol.GPMV, aes(Age, Cholesterol, color = Age))
d+
  geom_quasirandom(width = 0.2, size = 3, alpha = 0.8) +
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom = "errorbar", width=0.05)+
  stat_summary(fun=mean, geom = "point", size=15, shape=95) +
  scale_color_manual(values = pal4) +
  theme_classic() +
  labs(title = "Cholesterol on GPMV", y = "Concentration (uM)", x="") +
  stat_compare_means(aes(group = Age),label = "p.format") +
  theme(axis.text = element_text(colour = "black")) +
  theme(legend.position = "none") 
ggsave(filename = paste0("./Figure_Panels/EDFig.8e.pdf"), width = 3, height = 4, useDingbats = FALSE)

chol.med <- chol.GPMV %>% 
  group_by(Age) %>% 
  summarise(., MedianChol = median(Cholesterol))

chol.med
