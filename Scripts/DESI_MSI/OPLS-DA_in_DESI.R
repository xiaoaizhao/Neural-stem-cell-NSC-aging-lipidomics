
setwd(rstudioapi::getActiveProject())
rm(list=ls())

library(tidyverse)
library(ropls)

ssDT <- fread('./Output_Data/20210512_DESI_decomposition_single_sample.csv')

##====OPLS-DA between aNSCs and qNSCs from DESI====
all.DESI.Cell<- ssDT %>% 
  filter(!cell == "other") %>% 
  filter(!cell == "ki67") %>% 
  mutate(sample_cell = paste0(sampleID, "_", cell)) %>% 
  dplyr::select(peak, mean, sample_cell) %>% 
  pivot_wider(names_from = peak, values_from = mean) %>% 
  mutate(CellType = case_when(
    grepl("_gfap", sample_cell) ~ "qNSCs",
    grepl("_both", sample_cell) ~ "aNSCs",
  )) %>% 
  relocate(CellType, .after = "sample_cell")

cell.df <- all.DESI.Cell %>% 
  dplyr::select(CellType)

cell.mtx <- data.matrix(cell.df)

data.mtx <- all.DESI.Cell %>% 
  column_to_rownames(var = "sample_cell") %>% 
  dplyr::select(-CellType)

data.mtx <- data.matrix(data.mtx)
### ====Permutation test====
set.seed(12345)
Cell.oplsda <- opls(data.mtx, cell.mtx,
                        predI = 1, orthoI = 1, crossvalI = 6, permI = 100, scaleC = 'standard')
plot(Cell.oplsda, typeVc = "permutation") 

pR2Y <- Cell.oplsda@summaryDF$`pR2Y`
pQ2 <- Cell.oplsda@summaryDF$`pQ2`
cat("p-value for R²Y:", pR2Y, "\np-value for Q²:", pQ2)
# p-value for R²Y: 0.01 
# p-value for Q²: 0.01> 

perm_df <- as.data.frame(Cell.oplsda@suppLs$permMN)

ggplot(perm_df, aes(x = sim, y = `Q2(cum)`)) +
  geom_point(alpha = 0.65, size = 3) +
  geom_hline(yintercept = perm_df$`Q2(cum)`[perm_df$sim == 1], linetype = "dashed", color = "red") +
  labs(x = "Similarity", y = "Q²", title = "Permutation Test Results") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"))
ggsave(filename = "./Figure_Panels/EDFig.6d.pdf", width = 5, height = 5, useDingbats=FALSE)

### ====OPLS-DA====
o1 = as.data.frame(Cell.oplsda@orthoScoreMN)
p1 = as.data.frame(Cell.oplsda@scoreMN)

df.p <- bind_cols(o1, p1) %>% 
  rownames_to_column(var = "Sample") %>% 
  mutate(CellType = case_when(
    grepl("_gfap", Sample) ~ "qNSCs",
    # grepl("_ki67", Sample) ~ "NPCs",
    grepl("_both", Sample) ~ "aNSCs",
  )) %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

df.p$CellType <- factor(df.p$CellType, levels = c("aNSCs", "qNSCs"))
df.p$Age <- factor(df.p$Age, levels = c("Young", "Old"))

A <- ggplot(df.p, aes(p1, o1))
A + geom_point(aes( color = CellType), shape = 18, size = 4, alpha = 0.8)+
  scale_color_manual(values = c("#664692", "#63A66C"))+
  theme_classic() +
  theme(text=element_text(size = 13, face = "plain", colour = "black"))+
  theme(axis.text = element_text(colour = "black"))+
  ggtitle("DESI aNSC vs. qNSC - OPLS-DA")+
  xlab(paste0("Component 1: ", format(Cell.oplsda@modelDF["p1", "R2X"] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0("Orthogonal Component 1: ", format(Cell.oplsda@modelDF["o1", "R2X"] * 100,
                                  digits = 3), " % variance")) +
  theme(legend.position = "bottom")
ggsave(filename = "./Figure_Panels/EDFig.6c.pdf", width = 5, height = 5, useDingbats=FALSE)


##====OPLS-DA between young and old qNSCs from DESI====

Q.all.DESI<- ssDT %>% 
  filter(cell == "gfap") %>% 
  mutate(sample_cell = paste0(sampleID, "_", cell)) %>% 
  dplyr::select(peak, mean, sample_cell) %>% 
  pivot_wider(names_from = peak, values_from = mean) %>% 
  mutate(Age = case_when(
    grepl("^Y", sample_cell) ~ "Young", 
    grepl("^O", sample_cell) ~ "Old"
  )) %>% 
  mutate(CellType = case_when(
    grepl("_gfap", sample_cell) ~ "qNSCs",
  )) %>%
  relocate(c(Age, CellType), .after = "sample_cell")

var.mtx <- Q.all.DESI %>% 
  dplyr::select(Age)

var.mtx <- data.matrix(var.mtx)

data.mtx <- Q.all.DESI %>% 
  column_to_rownames(var = "sample_cell") %>% 
  dplyr::select(-c(Age, CellType))

Q.data.mtx <- data.matrix(data.mtx)

### ====Permutation test====
set.seed(1111)
Q.DESI.pls <- opls(Q.data.mtx, var.mtx, , predI = 1, orthoI = 1, crossvalI = 6, permI = 100, scaleC = 'standard')

plot(Q.DESI.pls, typeVc = "permutation") 

pR2Y <- Q.DESI.pls@summaryDF$`pR2Y`
pQ2 <- Q.DESI.pls@summaryDF$`pQ2`
cat("p-value for R²Y:", pR2Y, "\np-value for Q²:", pQ2)
# p-value for R²Y: 0.42 
# p-value for Q²: 0.98

perm_df <- as.data.frame(Q.DESI.pls@suppLs$permMN)

ggplot(perm_df, aes(x = sim, y = `Q2(cum)`)) +
  geom_point(alpha = 0.65, size = 3) +
  geom_hline(yintercept = perm_df$`Q2(cum)`[perm_df$sim == 1], linetype = "dashed", color = "red") +
  labs(x = "Similarity", y = "Q²", title = "Permutation Test Results") +
  theme_classic()
ggsave(filename = "./Figure_Panels/EDFig.6f.pdf", width = 5, height = 5, useDingbats=FALSE)

### ====OPLS-DA====
o1 = as.data.frame(Q.DESI.pls@orthoScoreMN)
p1 = as.data.frame(Q.DESI.pls@scoreMN)

df.p <- bind_cols(o1, p1) %>% 
  rownames_to_column(var = "Sample") %>% 
  mutate(Age = ifelse(grepl("^Y", Sample), "Young", "Old"))

df.p$Age <- factor(df.p$Age, levels = c("Young", "Old"))

pal4 <- c("cyan3", "magenta3")
A <- ggplot(df.p, aes(p1, o1))
A + geom_point(aes( color = Age), shape = 18, size = 4, alpha = 0.8)+
  scale_color_manual(values = pal4)+
  # scale_shape_manual(values = c(16, 17, 10)) +
  theme_classic() +
  theme(text=element_text(size = 13, face = "plain", colour = "black"))+
  theme(axis.text = element_text(colour = "black"))+
  ggtitle("DESI old vs. young qNSC - OPLS-DA")+
  xlab(paste0("Component 1: ", format(Q.DESI.pls@modelDF["p1", "R2X"] * 100,
                                  digits = 3), " % variance"))+
  ylab(paste0("Orthogonal Component 1: ", format(Q.DESI.pls@modelDF["o1", "R2X"] * 100,
                                  digits = 3), " % variance")) +
  theme(legend.position = "bottom")
ggsave(filename = "./Figure_Panels/EDFig.6e.pdf", width = 5, height = 5, useDingbats=FALSE)
