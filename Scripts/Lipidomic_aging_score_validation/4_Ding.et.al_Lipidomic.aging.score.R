
setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)

## ---------------------------------------------------------------------------------------------------------------------

load("./Output_Data/Ding.et.all.Lipid.zscore.on.aging.signature.Rdata")
load("./Output_Data/Ding.et.all.DB.zscore.on.aging.signature.Rdata")

z.avg.all <- full_join(Z.sc.sgnt.DB, Z.sc.sgnt.lpd, by = c("Samples", "Age", "Tissue", "Sex")) %>% 
  mutate(Mean.Z.signature = mean(c(Mean.Lpd.Z,  Mean.DB_Pct.Z)))

z.avg.all$Age <- factor(z.avg.all$Age, levels = c("3 weeks", "16 weeks", "59 weeks", "92 weeks"))

#' 
#' Plotting, add trend line between time points
## ---------------------------------------------------------------------------------------------------------------------
s.cor <- c("#fa715f", "#5b55d5")
z.avg.all.time <- z.avg.all %>% 
  mutate(TimePoint = as.numeric(substr(Age, 1, str_locate(Age, " ")-1)))

a <- ggplot(z.avg.all.time, aes(TimePoint, Mean.Z.signature, color = Sex))
a + geom_quasirandom(width = 0.3, alpha = 0.7) +
  geom_smooth(alpha = 0.5) +
  facet_wrap(~Tissue, nrow = 2) +
  theme_classic() +
  scale_color_manual(values = s.cor) +
  labs(title = "Lipidomic aging score (Ding et al.)", x = "Age (weeks)", y = "Mean Z score") +
  theme(axis.text = element_text(colour = "black"))
ggsave(filename = "./Figure_Panels/EDFig.11e.pdf", width = 8, height = 6, useDingbats = FALSE)

