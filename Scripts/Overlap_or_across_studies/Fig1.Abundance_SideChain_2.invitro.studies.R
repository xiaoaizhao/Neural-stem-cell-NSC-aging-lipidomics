## Need to run the following  script first
### 'Fig1.SideChain_2.invitro.studies.R'
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(ggrepel)

load("./Output_Data/SC.2outof2LC.invitro.Padj.Sig.features.Rdata")

load("./Output_Data/Invitro.SC.analysis.Rdata") 

Invitro.SC <- Invitro.SumSC %>% 
  mutate(Cla_SC = paste0(Class, "(", Sep.SC, ")")) %>% 
  mutate(Age = ifelse(grepl("^O", Samples), "Old", "Young")) %>% 
  group_by(Cla_SC) %>% 
  summarise(meanSumSC = mean(SumSC)) %>% 
  mutate(Exp = "In vitro")

load("./Output_Data/E2.SC.analysis.Rdata")

Exp2.ctrl.SC <- E2.SumSC %>% 
  mutate(KO = substr(Samples, nchar(Samples)-1, nchar(Samples))) %>% 
  filter(KO == "_N") %>% 
  mutate(Cla_SC = paste0(Class, "(", Sep.SC, ")")) %>% 
  mutate(Age = ifelse(grepl("^O", Samples), "Old", "Young")) %>% 
  group_by(Cla_SC) %>% 
  summarise(meanSumSC = mean(SumSC)) %>% 
  mutate(Exp = "In vitro #2")

### ==== 39 side chain features are too many to label, add text to only top 30% that go up and down with age ====
### ==== color code all significant SC that are identified in both in vitro studies ====
top30.df.hi <- padj.sig.db. %>%
  group_by(Cla_SC) %>%
  summarise(MeanES = unique(MeanES)) %>%
  mutate(percentile = percent_rank(MeanES)) %>%
  filter(percentile > 0.7) #18 side chain
top30.df.lo <- padj.sig.db. %>%
  filter(MeanES < 0) %>%
  group_by(Cla_SC) %>%
  summarise(MeanES = unique(MeanES)) %>%
  mutate(percentile = percent_rank(MeanES)) %>%
  filter(percentile < 0.3) #3 lipids
top30.ls <- c(top30.df.hi$Cla_SC, top30.df.lo$Cla_SC)

LC2 <- bind_rows(Invitro.SC, Exp2.ctrl.SC) %>% 
  rowwise() %>% 
  mutate(Sig.lbl = ifelse(Cla_SC %in% unique(padj.sig.db.$Cla_SC), "Significant", "Not significant")) %>% 
  mutate(Sig.txt = ifelse(Cla_SC %in% top30.ls, Cla_SC, "")) %>% 
  mutate(Log2.SCSum = log2(meanSumSC)) 

LC2$Sig.lbl <- factor(LC2$Sig.lbl, levels = c("Not significant" , "Significant"))
LC2$Exp <- factor(LC2$Exp, levels = c( "In vitro", "In vitro #2"))

LC2 <- LC2 %>% 
  arrange(Sig.lbl)

pal <- c("grey80", "orchid4")

a <- ggplot(LC2, aes(reorder(Cla_SC, Log2.SCSum), Log2.SCSum)) 
a+
  geom_point(aes(color = Sig.lbl), alpha = 0.85, size = 3)+
  # geom_point(aes(color = Sig.lbl), alpha = 0.75, size = 3)+
  theme_classic()+
  geom_text_repel(aes(label = Sig.txt), fontface = 'plain',
                  size = 3,colour = "black",
                  box.padding = unit(0.4, "lines"),
                  seed = 1234,
                  min.segment.length = 0,
                  max.overlaps = 55)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  labs(title = "Significant 2 out of 2 in vitro side chain abundance")+
  scale_x_discrete(expand = c(0.1,0.1)) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Log2 side chain concentration (uM) in respective class") +
  facet_wrap(~Exp,  nrow = 1)

ggsave(filename = "./Figure_Panels/EDFig.4c.pdf", width = 6, height = 5, useDingbats = FALSE)

t <- LC2 %>% 
  filter(Sig.lbl == "Significant") %>% 
  group_by(Exp) %>% 
  tally()
