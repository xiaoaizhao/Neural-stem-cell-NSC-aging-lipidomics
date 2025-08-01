##Prerequisite: Need to run `Lipids_invitro.fomatting.for.effect.size.correlation.R` first

## Correlation analysis on lipid effect size with age, on all datasets on Primary NSC cultures from LC-MS/MS and Lipidyzer lipidomics
## Correlation is done between 2 datasets at a time.


setwd(rstudioapi::getActiveProject())
rm(list = ls())
library(tidyverse)

load("./Output_Data/Ldz_Ef_Size_Age_qNSC.Rdata")# Lipidyzer
load("./Output_Data/Exp1_EfSz.age.ID_string.reFormat.dupTG.rmved.Rdata") # Primary culture #1
load("./Output_Data/Exp2.ctrl_EfSz.age.ID_string.reFormat.dupTG.rmved.Rdata") # Primary culture #2, control samples
load("./Output_Data/Exp3_EfSz.age.ID_string.reFormat.dupTG.rmved.Rdata")# Primary culture #3

## ==== Organize Lipidyzer data====
## Remove duplicated TGs before overlapping with other dataset, unique lipids from Lipidyzer is identified in script `Lipidyzer_data_remove_duplicated_lipids_before_overlapping.R`
load("./Output_Data/Lipidyzer_qNSC_for.LC-MS.ovlp.Rdata")
Ldz.ovlp <- Lipidyzer.Age.es.g %>% 
  filter(Lipid %in% Ldz.Qui.nodup$LipidIon) %>% 
  rowwise() %>% 
  mutate(ID_string = Ldz.Qui.nodup$ID_string[Ldz.Qui.nodup$LipidIon == Lipid])

conc.lpd.all4<- list(E3.fmt, 
                   E2.fmt, 
                   E1.fmt,
                   Ldz.ovlp) %>% 
  reduce(full_join, by = c("ID_string")) %>% #58 lipids
  rename("ES.Lipid_Exp3" = "es_g.x") %>% 
  rename("ES.Lipid_Exp2" = "es_g.y") %>% 
  rename("ES.Lipid_Exp1" = "es_g.x.x") %>% 
  rename("ES.Lipid_Lipidyzer" = "es_g.y.y") 

EfSz.df.all <- conc.lpd.all4 %>% 
  select(matches("ID_string|ES")) %>% 
  column_to_rownames(var = "ID_string")

## ====Plotting function====
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

## ====Make plot====
pdf(file = paste0("./Figure_Panels/EDFig.1j.pdf"))
pairs(EfSz.df.all,
      upper.panel = panel.cor_pval,
      lower.panel = panel.scatter_lm)
dev.off()
