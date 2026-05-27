
### Quantification of cholesterol concentration from GPMV lipidomics

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggbeeswarm)
library(ggpubr)

load("./Output_Data/GPMV.MsD.ClassSum.Rdata")
sum.df <- GPMV.class.sum %>%
  pivot_wider(., names_from = "Class", values_from = "ClassSum") %>%
  mutate(., POI = (PE+SM)/PC) %>%
  mutate(., Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

sum.df$Age <- factor(sum.df$Age, levels = c("Young", "Old"))

####====Cholesterol concentration====
chol.GPMV <- sum.df %>% 
  select(Sample, Cholesterol, Age)

chol.GPMV$Age <- factor(chol.GPMV$Age, levels = c("Young", "Old"))

pal4 <- c("cyan3", "magenta3")

d <- ggplot(chol.GPMV, aes(Age, Cholesterol, color = Age))
d+
  geom_quasirandom(width = 0.2, size = 3, alpha = 0.8) +
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom = "errorbar", width=0.05)+
  stat_summary(fun=median, geom = "point", size=15, shape=95) +
  scale_color_manual(values = pal4) +
  theme_classic() +
  labs(title = "Cholesterol on GPMV", y = "Concentration (uM)", x="") +
  stat_compare_means(aes(group = Age),label = "p.format") +
  theme(axis.text = element_text(colour = "black")) +
  theme(legend.position = "none") +
  ylim(0, max(chol.GPMV$Cholesterol+20))#+
  # geom_text_repel(aes(label = Sample))
ggsave(filename = paste0("./Figure_Panels/fig.S9H.pdf"), width = 3, height = 4, useDingbats = FALSE)

chol.med <- chol.GPMV %>% 
  group_by(Age) %>% 
  summarise(., MedianChol = median(Cholesterol))

chol.med
