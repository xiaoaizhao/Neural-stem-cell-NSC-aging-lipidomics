## Meta-analysis on double bond composition changes on 4 experiments
## Primary Culture #1 , Primary Culture #2, In vivo qNSC and GPMV
## Get DB features with significant summary effect size across all studies
##------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(data.table)
library(ggthemes)
library(rmeta)
library(RColorBrewer)

setwd(rstudioapi::getActiveProject())

## Load in effect size from each individual study, obtained from previous scripts
## ---------------------------------------------------------------------------------------------------------------
load('Output_Data/Ef_Size_DB_pct_InVitro.Rdata')
load('Output_Data/Ef_Size_DB_PCT_Exp2_all_KO.Rdata')
load('Output_Data/Ef_Size_DB_pct_Invivo.Rdata')
load('Output_Data/Ef_Size_DB_pct_GPMV.Rdata')


lipid_config_allExp_DT <- rbind(
  as.data.table(Invivo.DB.es.g)[,dataset:='In vivo sorted'],
  as.data.table(GPMV.DB.es.g)[,dataset:='GPMV'],
  as.data.table(LC.Invitro.DB.es.g)[,dataset:='Primary Culture'],
  as.data.table(Exp2.DB.es.g.allKO)[KO=="N"][,-13][,dataset:='Primary Culture #2']
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

summaryCastDT.DB <- dcast.data.table(
  lipid_config_allExp_DT[
    ,
    metaRunnerDT(es_g,
                 se_g)
    ,
    Cla_DB],Cla_DB~names,value.var='res')

summaryCastDT.DB[,zscore:=summary/se_summary]
summaryCastDT.DB[,pvalue:=2*pnorm(abs(zscore),lower.tail = F)]
save(summaryCastDT.DB, file = paste0("./Output_Data/Meta_DoubleBond_4_studies.Rdata"))


## ---------------------------------------------------------------------------------------------------------------
load("./Output_Data/Meta_DoubleBond_4_studies.Rdata")
higherOldSign <- summaryCastDT.DB[summary>0][order(pvalue)][pv_het>.05][n>2][pvalue<.05]

higherYngSign <- summaryCastDT.DB[summary<0][order(pvalue)][pv_het>.05][n>2][pvalue<.05]
DB.hi.old <- higherOldSign$Cla_DB
DB.hi.yng <- higherYngSign$Cla_DB

save(DB.hi.old, file = "./Output_Data/Meta_DB_signature_hi_in_Old_4_studies.Rdata")
save(DB.hi.yng, file = "./Output_Data/Meta_DB_signature_hi_in_Yng_4_studies.Rdata")


