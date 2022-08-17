library(boot)
library(csSAM)
library(parallel)
library(data.table)

csSAMcsFit <- function(inputPeakDT,inputCellDT){
  
  #convert DTs into matrices
  peakMat <- as.matrix(inputPeakDT[,-1])
  row.names(peakMat) <- inputPeakDT$sample
  
  cellMat <- as.matrix(inputCellDT[,-1])
  row.names(cellMat) <- inputCellDT$sample
  
  peakMat <- peakMat[order(row.names(peakMat)),]
  cellMat <- cellMat[order(row.names(cellMat)),]
  
  if(all(row.names(peakMat)!=row.names(cellMat))){warning("BUG!")}
  
  #run least square fit to get coefficients 
  #for cell-type specific values
  csfitOut <- csfit(cellMat,peakMat)
  
  #~~~~~~~
  expM <- csfitOut$ghat
  colnames(expM) <- colnames(peakMat)
  
  expDT <- melt(
    as.data.table(expM,
                  keep.rownames = T),
    id.vars = 'rn')
  
  #~~~~~~~
  serM <- csfitOut$se
  colnames(serM) <- colnames(peakMat)
  
  serDT <- melt(
    as.data.table(serM,
                  keep.rownames = T),
    id.vars = 'rn')
  
  #merge
  finalDT <- merge.data.table(expDT,
                              serDT,
                              by=c('rn','variable'))
  setnames(finalDT,c('cell','peak','mean','se'))
  finalDT[,n:=nrow(peakMat)]
  
  return(finalDT)
}

computeT <- function(m1,m2,se1,se2,n1,n2){
  
  s1  <- se1*(sqrt(n1))
  s2  <- se2*(sqrt(n2))
  
  tn  <- m1-m2
  sp  <- sqrt(((n1-1)*(s1**2)+(n2-1)*(s2**2))/(n1+n2-2))
  td  <- sp*sqrt(1/n1 + 1/n2)
  
  return(c(tn/td,n1+n2-2))
}

computeWT <- function(m1,m2,se1,se2,n1,n2){
  
  s1  <- se1*(sqrt(n1))
  s2  <- se2*(sqrt(n2))
  
  tn  <- m1-m2
  td  <- sqrt((s1**2/n1) + (s2**2/n2))
  
  dof_n <- (((s1**2)/n1) + ((s2**2)/n2))**2
  dof_d <- (s1**4)/((n1**2)*(n1-1)) + (s2**4)/((n2**2)*(n2-1))
  
  return(c(tn/td,round(dof_n/dof_d,0)))
}

computePfromT <- function(t,dof){
  2*pt(abs(t),dof,lower.tail = F)
}


