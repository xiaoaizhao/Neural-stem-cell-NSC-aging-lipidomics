## Correlation analysis on lipid effect size with age, between in vivo isolated qNSC and 3 Primary NSC cultures from LC-MS/MS
## Correlation is done between 2 datasets at a time.

setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)
source("./Scripts/Function_scripts/Effect_size_functions.R")
load("./Output_Data/Ef_Size_Lipid_Age_qNSC_InVitro_LC.Rdata") # Primary culture #1
load("./Output_Data/Ef_Size_CONC.Lipid_Exp2_all_KO.Rdata") # Primary culture #2
load("./Output_Data/Exp3_Qui_Age_ES.Rdata") # Primary culture #3
load("./Output_Data/Ef_Size_CONC.Lipid_InVivo.Rdata") # In vivo isolated qNSCs


E1.n <- ether.rename(InVitro.lpd.es.g)
E2.n <- ether.rename(Exp2.CONC.lpd.es.g.allKO) %>% 
  filter(KO == "N")
E3.n <- E3.Q.AG.ES
invivo.n <- ether.rename(Invivo.CONC.lpd.es.g)

### All 4 experiments together 1X1

conc.lpd.all4<- list(E1.n, 
                     E2.n, 
                     E3.n,
                     invivo.n) %>% 
  reduce(full_join, by = c("LipidIon")) %>% #58 lipids
  rename("ES.Lipid_Exp1" = "es_g.x") %>% 
  rename("ES.Lipid_Exp2" = "es_g.y") %>% 
  rename("ES.Lipid_Exp3" = "es_g.x.x") %>% 
  rename("ES.Lipid_Invivo" = "es_g.y.y") 

EfSz.df.all <- conc.lpd.all4 %>% 
  select(matches("LipidIon|ES")) %>% 
  column_to_rownames(var = "LipidIon")


# Custom upper panel (correlation + p-value)
panel.cor_pval <- function(x, y, digits = 2, cex.cor = 1.3, ...) {
  test <- cor.test(x, y, method = "pearson", na.action = na.omit)
  cor_val <- round(test$estimate, digits)
  p_val <- ifelse(test$p.value < 0.01, "<0.01", round(test$p.value, digits))
  
  # Line 1: Correlation coefficient (r)
  text(mean(x, na.rm = TRUE), mean(y, na.rm = TRUE) + 0.1 * diff(range(y, na.rm = TRUE)),  # Shift upward for 1st line
       bquote(italic(r) == .(cor_val)), 
       cex = cex.cor)
  
  # Line 2: P-value
  text(mean(x, na.rm = TRUE), mean(y, na.rm = TRUE) - 0.1 * diff(range(y, na.rm = TRUE)),  # Shift downward for 2nd line
       bquote(italic(p) == .(p_val)), 
       cex = cex.cor,
       col = ifelse(test$p.value < 0.05, "red", "black"))
}

# Custom lower panel (scatterplot + regression line)
mycolor <- adjustcolor("grey25", alpha.f = 0.6)
panel.scatter_lm <- function(x, y, col = mycolor, pch = 19, cex = 1, 
                             col.smooth = "blue3", lty.smooth = 1, lwd = 1.5, ...) {
  points(x, y, col = col, pch = pch, cex = cex)
  abline(lm(y ~ x), col = col.smooth, lty = lty.smooth, lwd = lwd)
}

# Plot
pdf(file = paste0("./Figure_Panels/EDFig.5f.pdf"))
pairs(EfSz.df.all,
      upper.panel = panel.cor_pval,
      lower.panel = panel.scatter_lm)
dev.off()
