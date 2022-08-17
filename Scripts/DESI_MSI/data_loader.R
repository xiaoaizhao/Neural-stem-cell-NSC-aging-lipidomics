rm(list=ls())
library(data.table)

peakLoader <- function(peakFolder = './Input_Data/Cleaned_title_all_spectra'){
  
  #aggregate all files
  massDT <- rbindlist(
    lapply(
      list.files(peakFolder,full.names = T),
      function(i)
        fread(i,skip = 1)[,name:=gsub(".txt","",gsub(".*/","",i))]))
  setnames(massDT,c('row','raw_peak','intensity','sample'))
  
  #organize data by subject and slice
  massDT[,  slice:=gsub(".*_","",sample)]
  massDT[,subject:=gsub("_(1|2|3|4|5|6|7|8|9|0)+$","",sample)]
  
  return(massDT)
}

propLoader <- function(propCfile = './Input_Data/CellProportions/Censored_allout_Xiaoai_percentages_9.11.18.csv',
                       propUfile = './Input_Data/CellProportions/Uncensored_allout_Xiaoai_percentages_9.11.18.csv'){
  
  propsCensored   <- fread(propCfile)
  propsUncensored <- fread(propUfile)
  
  propsUncensored[,highQ:=0]
  propsUncensored[file %in% propsCensored$file,highQ:=1]
  propsUncensored[,sample:=gsub(".tif","",file)]
  
  return(propsUncensored)
}

peakNormalizer <- function(inputPeakDT,
                           intStIntF = './Input_Data/CellProportions/InternalStandard_intensity_allSample_091218.csv'){
  return(merge.data.table(inputPeakDT,
                          fread(intStIntF),
                          by.x='sample',
                          by.y='Image')[,intensity:=intensity/`IntStd(max)`])
}

roundingMassPeak <- function(inputPeakDT){
  
  #perform rounding of mass peak features
  inputPeakDT[,  roundedP:=round(raw_peak)]
  inputPeakDT[,roundedP_1:=round(raw_peak,1)]
  inputPeakDT[,roundedP_2:=round(raw_peak,2)]
  
  return(inputPeakDT)
}

