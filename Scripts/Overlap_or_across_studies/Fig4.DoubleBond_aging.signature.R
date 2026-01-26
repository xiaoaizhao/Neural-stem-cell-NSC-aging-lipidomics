### Meta-analysis on double bond composition changes to obtain lipidomic aging signature
##------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(data.table)
library(rmeta)

setwd(rstudioapi::getActiveProject())

## Load in effect size from each individual study, obtained from previous scripts
## ---------------------------------------------------------------------------------------------------------------
load('Output_Data/DBPct.MsD.Invitro_Qui_Age_ES.Rdata')
load('Output_Data/DB.MsD.Age.ES_Exp2_all_KO.Rdata')
load('Output_Data/DBPct.MsD.Invivo_Age_ES.Rdata')
load('Output_Data/DBPct.MsD.GPMV_Age_ES.Rdata')


lipid_config_allExp_DT <- rbind(
  as.data.table(Invivo.DB.AG.ES)[,dataset:='In vivo'],
  as.data.table(GPMV.DB.AG.ES)[,dataset:='GPMV'],
  as.data.table(Invitro.Qui.MsD.CONC.DB.es.g)[,dataset:='In vitro'],
  as.data.table(Exp2.MsD.DB.es.g.allKO)[KO=="N"][,-13][,dataset:='In vitro Experiment #2']
) #194


## ---------------------------------------------------------------------------------------------------------------
lipid_config_allExp_DT <- lipid_config_allExp_DT[!is.na(es_g)] #178


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
summaryCastDT.DB <- summaryCastDT.DB %>% 
  mutate(padj = p.adjust(pvalue, method = "fdr"))
DB.summary <- summaryCastDT.DB
save(DB.summary, file = "./Output_Data/DB.summary.from.meta.analysis.Rdata")
## ---------------------------------------------------------------------------------------------------------------
higherOldSign <- summaryCastDT.DB[summary>0][order(pvalue)][pv_het>.05][n>1][pvalue<.15] 


conc.DB.hi.old <- higherOldSign$Cla_DB
conc.DB.hi.old <- conc.DB.hi.old[!grepl("PI:4", conc.DB.hi.old)] #5 features

save(conc.DB.hi.old, file = "./Output_Data/Meta_DB_signature.Rdata")


