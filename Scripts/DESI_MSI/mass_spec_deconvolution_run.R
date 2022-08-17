################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#DESI in-vivo Lipid Mass Spec Decomposition
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#R functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls())
library(data.table)
source('./Scripts/DESI_MSI/data_loader.R')
source('./Scripts/DESI_MSI/deconvolution_functions.R')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#R Main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Load input data [with peak normalization]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#aggregate all files
massDT          <- roundingMassPeak(peakNormalizer(peakLoader()))
propsUncensored <- propLoader()

#merge
mass_props_DT <- merge.data.table(massDT,propsUncensored,by='sample')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Combine data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

peakDT <- dcast.data.table(mass_props_DT,
                           sample~roundedP_1,
                           value.var = 'intensity',
                           fun.aggregate = sum)
peakMat <- as.matrix(peakDT[,-1])
row.names(peakMat) <- peakDT$sample

cellDT <- unique(mass_props_DT[,c('sample',
                                  'percent_ki67_alone',
                                  'percent_gfap_alone',
                                  'percent_both')])
setnames(cellDT,c('sample','ki67','gfap','both'))
cellDT[,other:=1-(ki67+gfap+both)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#csSAM: single subject level (mouse)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#get all subject IDs [which map to 10 samples]
subjectsIDs <- unique(
  sapply(peakDT[,.N,sample]$sample,
         function(i)
           paste(strsplit(i,"_")[[1]][1:2],collapse =  "_")))

singleSubjectDecompositionDT <- lapply(
  subjectsIDs,
  function(sampleID){
    
    #run csSAM
    csSAMfitTemp <- csSAMcsFit(
      peakDT[grep(sampleID,sample)],
      cellDT[grep(sampleID,sample)])
    
    #label them
    csSAMfitTemp[,sampleID:=sampleID]
    
    return(
      csSAMfitTemp
    )
  })

singleSubjectDecompositionDT <- rbindlist(singleSubjectDecompositionDT)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#compute Effect Size
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

singleSubjectDecompositionDT[,class:=1]
singleSubjectDecompositionDT[grep("Y",sampleID),class:=0]

#compute effect size [es]
#old vs young comparison 
#(pos ES -> higher in old)
#(neg ES -> higher in young)
#adapted from MetaIntegrator (CRAN)
es.g <- function(valueV,classV){
  
  diff <- mean(valueV[classV==1]) - mean(valueV[classV==0])
  
  sd1  <- sd(valueV[classV==0])
  sd2  <- sd(valueV[classV==1])
  n1   <- length(valueV[classV==0])
  n2   <- length(valueV[classV==1])
  
  sp   <- sqrt( ( (n1-1)*sd1^2 + (n2-1)*sd2^2 )/( n1 + n2 - 2 ) )
  cf   <- 1 - 3/( 4*(n1 + n2) - 9 ) # small sample size correction
  
  g    <- cf * diff/sp
  se.g <- sqrt( (n1+n2)/(n1*n2) + 0.5*g^2 /(n1+n2-3.94))
  
  return(c(g,se.g))
}

#format data
esDT <- singleSubjectDecompositionDT[
  ,es.g(mean,
        class),
  .(peak,cell)]
esDT[,names:=rep(c('es','se'),.N/2)]
esDT <- dcast.data.table(esDT,peak+cell~names,value.var = 'V1')

#compute p-value from effect size
esDT[,pv:=2*pnorm(abs(es/se),lower.tail = F)]
esDT <- esDT[order(pv)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Save output file into a text
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fwrite(file = './Output_Data/20210512_DESI_decomposition_single_sample.csv',
       singleSubjectDecompositionDT,
       quote = F)

fwrite(file = './Output_Data/20210512_DESI_decomposition_ES_oldVsYoung.csv',
       esDT,
       quote = F)
