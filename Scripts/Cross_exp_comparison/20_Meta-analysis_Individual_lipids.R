## Meta-analysis on individual lipids on 4 experiments
## Primary Culture #1 , Primary Culture #2, In vivo qNSC and GPMV
## Get lipids with significant summary effect size across all studies
##------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(data.table)
library(ggthemes)
library(rmeta)
library(RColorBrewer)
library(rstatix)

setwd(rstudioapi::getActiveProject())

## Load in effect size from each individual study, obtained from previous scripts
## ---------------------------------------------------------------------------------------------------------------
load("./Output_Data/Ef_Size_Lipid_Age_qNSC_InVitro_LC.Rdata")
load("./Output_Data/Ef_Size_Lipid_Exp2_all_KO.Rdata")
load("./Output_Data/Ef_Size_Lipid_GPMV.Rdata")
load("./Output_Data/Ef_Size_Lipid_InVivo.Rdata")


lipid_config_allExp_DT <- rbind(
  as.data.table(Invivo.lpd.es.g)[,dataset:='In vivo sorted'],
  as.data.table(GPMV.lpd.es.g)[,dataset:='GPMV'],
  as.data.table(InVitro.lpd.es.g)[,dataset:='Primary Culture'],
  as.data.table(Exp2.lpd.es.g.allKO)[KO=="N"][,-13][,dataset:='Primary Culture #2']
)


## ---------------------------------------------------------------------------------------------------------------
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

summaryCastDT <- dcast.data.table(
  lipid_config_allExp_DT[
    ,
    metaRunnerDT(es_g,
                 se_g)
    ,
    LipidIon],LipidIon~names,value.var='res')

summaryCastDT[,zscore:=summary/se_summary]
summaryCastDT[,pvalue:=2*pnorm(abs(zscore),lower.tail = F)]
save(summaryCastDT, file = paste0("./Output_Data/Meta_lipids_4_studies.Rdata"))


## ---------------------------------------------------------------------------------------------------------------
load("./Output_Data/Meta_lipids_4_studies.Rdata")
higherOldSign <- summaryCastDT[summary>0][order(pvalue)][pv_het>.05][n>2][pvalue<.05] #22 lipids
higherYngSign <- summaryCastDT[summary<0][order(pvalue)][pv_het>.05][n>2][pvalue<.05] #20 lipids

lpd.hi.old <- higherOldSign$LipidIon
lpd.hi.yng <- higherYngSign$LipidIon

save(lpd.hi.old, file = "./Output_Data/Meta_Lipid_signature_hi_in_Old_4_studies.Rdata")
save(lpd.hi.yng, file = "./Output_Data/Meta_Lipid_signature_hi_in_Yng_4_studies.Rdata")


