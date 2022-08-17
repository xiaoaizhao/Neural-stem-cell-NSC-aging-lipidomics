## Objective: assess age-related changes in lipids that are detected by LC-MS datasets vs. DESI-MSI
## 3 LC-MS datasets are used in total, Primary Culture #1, Primary Culture #2 and In vivo qNSC lipidomic dataset
## Steps:
# 1. Perform meta-analysis on all 3 LC-MS datasets
# 2. Overlap lipids that are detected and validated on DESI-MSI
# 3. Perform correlation analysis
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(data.table)
library(ggthemes)
library(rmeta)
library(mgcv)
library(RColorBrewer)
library(ggrepel)

setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load("./Output_Data/Ef_Size_Lipid_Age_qNSC_InVitro_LC.Rdata")
load("./Output_Data/Ef_Size_Lipid_Exp2_all_KO.Rdata")
load("./Output_Data/Ef_Size_Lipid_InVivo.Rdata")


lipid_config_allExp_DT <- rbind(
  as.data.table(Invivo.lpd.es.g)[,dataset:='In vivo sorted'],
  as.data.table(InVitro.lpd.es.g)[,dataset:='Primary Culture'],
  as.data.table(Exp2.lpd.es.g.allKO)[KO=="N"][,-13][,dataset:='Primary Culture #2']
)

## Meta-analysis on 3 datasets
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
lipid_config_allExp_DT <- lipid_config_allExp_DT[!is.na(es_g)]


metaRunnerDT <- function(esV,
                         seV){
  
  metaR <- meta.summaries(
    esV,
    seV,
    method = 'random')
  
  outV <- c(metaR$summary,
            metaR$se.summary,
            metaR$het[3],
            length(metaR$effects))
  return(
    data.table(
      res=outV,
      names=c('summary',
              'se_summary',
              'pv_het',
              'n')))
}

summaryCastDT.noGPMV <- dcast.data.table(
  lipid_config_allExp_DT[
    ,
    metaRunnerDT(es_g,
                 se_g)
    ,
    LipidIon],LipidIon~names,value.var='res')

summaryCastDT.noGPMV[,zscore:=summary/se_summary]
summaryCastDT.noGPMV[,pvalue:=2*pnorm(abs(zscore),lower.tail = F)]
save(summaryCastDT.noGPMV, file = "./Output_Data/Meta_lipids_3_studies_no_GPMV.Rdata")


## Add m/z value to lipids after meta-analysis. If a lipid is detected in multiple studies, will use the mean m/z and convert all ion to [M-H]- to overlap with DESI-MSI data
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load("./Output_Data/Meta_lipids_3_studies_no_GPMV.Rdata")

ef_g_noGPMV <- summaryCastDT.noGPMV %>%
  filter(., !is.na(summary))

##match primary culture #1
load("./Output_Data/2_Ion+FA_clean_373_lipids.Rdata") ##lipid data frame with m/z for each lipid species

Sum_ef_g.mz <-  add.mz(FA_unique373, m.z, LipidIon, ef_g_noGPMV) %>%
  rename(., "m.z_Exp1" = "m.z") 

##match primary culture #2
load("./Output_Data/Exp2_Dup_ID_rmv+dup_isomer_rmv_694_lipids.Rdata")
Sum_ef_g.mz <- Sum_ef_g.mz %>%
  rename("LipidIon" = "LipidIonFix") 

Sum_ef_g.mz2 <- add.mz(df.rmv.isomer, CalcMz, LipidIon, Sum_ef_g.mz) %>%
  rename(., "m.z_Exp2" = "CalcMz") 

##match in vivo sorted
load("./Output_Data/Invivo_Ion+FA_clean_130_lipids.Rdata")
Sum_ef_g.mz2 <- Sum_ef_g.mz2 %>%
  rename("LipidIon" = "LipidIonFix") 
Sum_ef_g.mz3 <- add.mz(FA_unique130, m.z, LipidIon, Sum_ef_g.mz2) %>%
  rename(., "m.z_invivo" = "m.z")
Sum_ef_g.mz3 <- Sum_ef_g.mz3 %>%
  group_by(LipidIonFix) %>%
  group_modify(~{
    .x %>%
      mutate(., MZ_avg = mean(c(m.z_Exp1, m.z_Exp2, m.z_invivo), na.rm = T))
  }) %>%
  rename("LipidIon" = "LipidIonFix") 
  

Sum_ef_g.NegH <- Conv.NegH(Sum_ef_g.mz3) %>%
  ungroup() %>%
  rename(., "Summary_Ef_size" = "summary") %>%
  select(c("peak", "Summary_Ef_size", "LipidIon"))

## Input list of lipids that were validated by a separate tandem MS for structural identification, only use this set of lipids for correlation analysis
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DESI.ls <- read.csv("./Input_Data/DESI_validated_lipids_list_072922.csv", stringsAsFactors = F)

# reformat to string for subsetting
d <- str_c(DESI.ls$Ion, collapse = ",")

DESI.ls.c <-DESI.ls %>% 
  mutate(., Ion = Lbr(Ion)) %>% 
  mutate(., Ion = Rbr(Ion))

ID.list <- str_c(DESI.ls.c$Ion, collapse = ",") 
ID.list <- str_replace_all(ID.list, ",", "|")


Sum.ef.DESI <- Sum_ef_g.NegH %>% 
  mutate(., NoIon = substr(LipidIon, 1, str_locate(LipidIon, "\\)"))) %>%
  filter(., str_detect(NoIon, ID.list)) 

## Read in effect size matrix of DESI dataset
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
DESI <- read.csv("./Output_Data/20210512_DESI_decomposition_ES_oldVsYoung.csv", stringsAsFactors = F)

DESI_GFAP <- DESI %>%
  filter(., cell == "gfap") %>%
  select(c("peak", "es"))

DESIv3Exps <- inner_join(DESI_GFAP, Sum.ef.DESI, by = "peak") 

## Plot for correlation
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

a <- ggplot(DESIv3Exps, aes(Summary_Ef_size, es))
a+geom_smooth(method = 'lm', alpha = 0.3, colour = "blue4", fill='lightskyblue',linetype=1)+
  geom_point(size = 3, colour = "grey25")+
  stat_cor(method = "pearson", label.x = -2.5, label.y = 0.5)+
  theme_classic()+
  labs(title = "DESI vs. all LC" , 
       x = "Summary effect size from 3 LC experiments", 
       y = "DESI effect size")+
  geom_text_repel(aes(label = NoIon), fontface = 'plain',
                  size = 3,colour = "black",
                  xlim = c(min(DESIv3Exps$Summary_Ef_size)-0.2, 
                           max(DESIv3Exps$Summary_Ef_size)+0.2),
                  box.padding = unit(0.55, "lines"),
                  seed = 1234,
                  min.segment.length = 0,
                  max.overlaps = 10)+
  theme(text = element_text(size = 11, face = "plain"))+
  theme(axis.text = element_text(size = 11, face = "plain", colour = "black"))+
  theme(legend.position = "none")
ggsave(filename = "./Figure_Panels/Fig_2e.pdf", width = 5, height = 5,
       useDingbats=FALSE)


