### Meta-analysis on lipid changes to obtain lipidomic aging signature
##------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(data.table)
library(rmeta)

setwd(rstudioapi::getActiveProject())
source("./Scripts/Function_scripts/Pre-processing_functions.R")
## Load in effect size from each individual study, obtained from previous scripts
## ---------------------------------------------------------------------------------------------------------------
load("./Output_Data/Lpd.MsD.Invitro_Qui_Age_ES.w.RT.Rdata")
load("./Output_Data/Lpd.MsD.Age.ES_Exp2_all_KO.w.RT.Rdata")
load("./Output_Data/Lpd.MsD.GPMV_Age_ES.Rdata")
load("./Output_Data/Lpd.MsD.Invivo_Age_ES.w.RT.Rdata")

GPMV.meta <- MsD.lpd.rmv.abc(GPMV.lpd.ES.w.RT) 

Invivo.meta <- MsD.lpd.rmv.abc(Invivo.lpd.ES.w.RT)

KO.meta <- MsD.lpd.rmv.abc(Exp2.lpd.es.g.MsD.allKO.wRT)

Invitro.meta <- MsD.lpd.rmv.abc(Invitro.Q.lpd.ES.w.RT)

lipid_config_allExp_DT <- rbind(
  as.data.table(Invivo.meta)[,dataset:='In vivo sorted'],
  as.data.table(GPMV.meta)[,dataset:='GPMV'],
  as.data.table(Invitro.meta)[,dataset:='In vitro'],
  as.data.table(KO.meta)[KO=="N"][,-14][,dataset:='In vitro Experiment #2']
)

## ---------------------------------------------------------------------------------------------------------------
lipid_config_allExp_DT <- lipid_config_allExp_DT[!is.na(es_g)] #1002


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
summaryCastDT <- summaryCastDT %>% 
  mutate(padj = p.adjust(pvalue, method = "fdr"))
Lipid.summary <- summaryCastDT
save(Lipid.summary, file = "./Output_Data/Lipid.summary.from.meta.analysis.Rdata")
## ---------------------------------------------------------------------------------------------------------------

higherOldSign <- summaryCastDT[summary>0][order(pvalue)][pv_het>.05][n>2][pvalue<0.1] #25 

conc.lpd.hi.old <- higherOldSign$LipidIon
save(conc.lpd.hi.old, file = "./Output_Data/Meta_Lipid_signature.Rdata")

