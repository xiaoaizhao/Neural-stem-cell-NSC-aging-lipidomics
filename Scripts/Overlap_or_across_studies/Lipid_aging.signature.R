## Meta-analysis on lipid changes on 4 experiments: Primary culture #2, Primary culture #3, in vivo and GPMV

## Get lipids with significant summary effect size across all studies

##------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(data.table)
library(rmeta)

setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Effect_size_functions.R")
## Load in effect size from each individual study, obtained from previous scripts
## ---------------------------------------------------------------------------------------------------------------
load("./Output_Data/Ef_Size_CONC.Lipid_Exp2_all_KO.Rdata")
load("./Output_Data/Exp3_Qui_Age_ES.conc.lipids.Rdata")
load("./Output_Data/Ef_Size_Conc.Lipid_GPMV.Rdata")
load("./Output_Data/Ef_Size_CONC.Lipid_InVivo.Rdata")

GPMV.meta <- ether.rename(GPMV.conc.lpd.es.g)
Invivo.meta <- ether.rename(Invivo.CONC.lpd.es.g)
KO.meta <- ether.rename(Exp2.CONC.lpd.es.g.allKO)
E3.meta <- E3.Q.conc.AG.ES


lipid_config_allExp_DT <- rbind(
  as.data.table(Invivo.meta)[,dataset:='In vivo sorted'],
  as.data.table(GPMV.meta)[,dataset:='GPMV'],
  as.data.table(E3.meta)[,dataset:='Primary culture #3'],
  as.data.table(KO.meta)[KO=="N"][,-13][,dataset:='Primary culture #2']
)

## ---------------------------------------------------------------------------------------------------------------
lipid_config_allExp_DT <- lipid_config_allExp_DT[!is.na(es_g)]


metaRunnerDT <- function(esV,
                         seV){
  # set.seed(1314)
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
Lipid.summary <- summaryCastDT
save(Lipid.summary, file = "./Output_Data/Lipid.summary.from.meta.analysis.Rdata")
## ---------------------------------------------------------------------------------------------------------------

higherOldSign <- summaryCastDT[summary>0][order(pvalue)][pv_het>.05][n>1][pvalue<.05] #26/28  conc/all lipids

conc.lpd.hi.old <- higherOldSign$LipidIon

save(conc.lpd.hi.old, file = "./Output_Data/Meta_Lipid_signature.Rdata")



