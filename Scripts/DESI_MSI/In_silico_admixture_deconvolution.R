## ----------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
library(ragg)
library(tidyverse)
library(data.table)
source('./Scripts/DESI_MSI/data_loader.R')
source('./Scripts/DESI_MSI/deconvolution_functions.R')


## ----------------------------------------------------------------------------------------------------------------------------------
load("./Output_data/LC_Invitro_for_admixture.Rdata")
Invitro <- LC.AQ %>% 
  select(-c(Class, SideChain))



## ----------------------------------------------------------------------------------------------------------------------------------
load("./Output_Data/Neuron_Oligo_for_admixture.Rdata")


## ----------------------------------------------------------------------------------------------------------------------------------
col.ls <- c("#7D65A5", "#81B486", "darkorange4", "deepskyblue4")


## ----------------------------------------------------------------------------------------------------------------------------------
invitroDT <- melt.data.table(as.data.table(Invitro),
                             id.vars = 'ID_string')
invitroDT[,status:='Quiescent']
invitroDT[grep("aNSC-A",variable),status:='Activated']
invitroDT[,subj:=gsub("_aNSC-A|_qNSC-Q","",variable)]


## ----------------------------------------------------------------------------------------------------------------------------------
NEUoligoDT <- melt.data.table(as.data.table(BR.Neu.Oligo),
                              id.vars = 'ID_string')
NEUoligoDT[,status:='Neuron']
NEUoligoDT[grep("Oligo",variable),status:='Oligodendrocytes']
NEUoligoDT[,subj:=gsub("Oligodendrocytes|Neurons","",variable)]


## ----------------------------------------------------------------------------------------------------------------------------------
vectorNospace <- function(inVec){
  sapply(inVec,function(i)
    gsub(" +","",i))
}

listSortCollapse <- function(inlist){
  paste(sort(vectorNospace(inlist)),collapse = ',')
}

sortCollapse <- function(instring){
  listSortCollapse(strsplit(instring,",")[[1]])
}


NEUoligoDT[,newID:=sapply(NEUoligoDT$ID_string,
                          function(i)
                            listSortCollapse(i))]


invitroDT[,newID:=sapply(invitroDT$ID_string,
                         function(i)
                           listSortCollapse(i))]


combDT <- rbind(
  NEUoligoDT,invitroDT)[,mean(value,na.rm=T),
                        .(variable,newID)][!is.na(V1)]



## ----------------------------------------------------------------------------------------------------------------------------------
#for every subject get all cells, then generate n admixtures
mixEngine <- function(inputCastDT){
  
  #generate random proportion
  propV <- sample(1:100,ncol(inputCastDT)-1,replace = T)
  propV <- propV/sum(propV)
  
  #compute mix
  return(data.table(
    rn=inputCastDT$rn,
    mix=as.matrix(inputCastDT[,-1]) %*% propV,
    Activated = propV[1],
    Neurons = propV[2],
    Oligodendrocytes = propV[3],
    Quiescent = propV[4]#this is hardcoded now,fix later
  ))
}

mixSbjRunner <- function(inputDT,
                         subject,
                         mixNumber){
  
  #cast data.table
  castDT <- dcast.data.table(inputDT[subj==subject],rn~status)
  
  #generate proportion vector
  return(
    rbindlist(lapply(1:mixNumber,function(i)
      mixEngine(castDT)[,n:=i])
    )[,sbj:=subject])
}

mixSbjWrapper <- function(inputDT,...){
  
  #for each subject compute admixtures
  return(
    rbindlist(
      lapply(inputDT[,.N,subj]$subj,
             function(sID){
               mixSbjRunner(inputDT,sID,...)
             })))
}


## ----------------------------------------------------------------------------------------------------------------------------------
admixture_input_DT <- copy(
  combDT[newID %in% combDT[,.N,newID][order(N)][N==37]$newID]
)

admixture_input_DT[,Lsubj:=gsub(
  "Activated|Quiescent|Neurons|Oligodendrocytes",
  "",
  variable)] #1776


## ----------------------------------------------------------------------------------------------------------------------------------
set.seed(3333)

invitroSamplesA <- sample(admixture_input_DT[grep("_aNSC-A",variable),
                                            .N,Lsubj]$Lsubj,9)

invitroSamplesQ <- sample(admixture_input_DT[grep("_qNSC-Q",variable),
                                            .N,Lsubj]$Lsubj,9)

neuronSamples  <- sample(admixture_input_DT[grep("Neuron",variable),
                                            .N,Lsubj]$Lsubj,9)

oligosSamples  <- sample(admixture_input_DT[grep("Oligo",variable),
                                            .N,Lsubj]$Lsubj,9)

match_id_DT <- as.data.table(rbind(
  cbind(invitroSamplesA,invitroSamplesA),
  cbind(invitroSamplesQ,invitroSamplesA),
  cbind(neuronSamples,invitroSamplesA),
  cbind(oligosSamples,invitroSamplesA))
)
setnames(match_id_DT,c("Lsubj","subj"))



## ----------------------------------------------------------------------------------------------------------------------------------
admixture_input_DT <- admixture_input_DT[match_id_DT,o=.(Lsubj)]


## ----------------------------------------------------------------------------------------------------------------------------------
admixture_input_DT[,status:='Quiescent']
admixture_input_DT[grep("aNSC-A",variable),status:='Activated']

admixture_input_DT[grep("Neur",variable),status:='Neurons']
admixture_input_DT[grep("Oligo",variable),status:='Oligodendrocytes']


## ----------------------------------------------------------------------------------------------------------------------------------
admixture_input_DT[,value:=scale(V1),.(newID,status)]
admixture_input_DT[,rn:=newID]


## ----------------------------------------------------------------------------------------------------------------------------------
set.seed(3333)
invitro_neurons_MixDT <- mixSbjWrapper(admixture_input_DT,mixNumber=10)


## ----------------------------------------------------------------------------------------------------------------------------------
mixPropDT <- unique(
  invitro_neurons_MixDT[,c('Activated',
                           'Neurons',
                           'Oligodendrocytes',
                           'Quiescent',
                           'n','sbj')])
mixPropDT[,sample:="AAA"]

for(i in 1:nrow(mixPropDT)){
  mixPropDT[i,sample:=paste(sbj,n,sep = "_")]
}

mixPropDT <- mixPropDT[,c('sample',
                          'Activated',
                          'Neurons',
                          'Oligodendrocytes',
                          'Quiescent')]


## ----------------------------------------------------------------------------------------------------------------------------------
mixPeakDT <- unique(invitro_neurons_MixDT[,c('rn','mix.V1','n','sbj')])
mixPeakDT[,sample:="AAA"]

for(i in 1:nrow(mixPeakDT)){
  mixPeakDT[i,sample:=paste(sbj,n,sep = "_")]
}
mixPeakDT <- dcast.data.table(mixPeakDT,sample~rn,value.var = 'mix.V1')


## ----------------------------------------------------------------------------------------------------------------------------------
csSAMfitYadMIX <- csSAMcsFit(mixPeakDT[grep("Y",sample)],
                             mixPropDT[grep("Y",sample)])

csSAMfitOadMIX <- csSAMcsFit(mixPeakDT[grep("O",sample)],
                             mixPropDT[grep("O",sample)])


## ----------------------------------------------------------------------------------------------------------------------------------
#functionalized single sample deconvolution runner

singleSampleDecRunner <- function(inputPeakDT,
                                  inputPropDT){
  
  subjV <- unique(
    sapply(inputPeakDT[,.N,sample]$sample,
           function(i)
             paste(strsplit(i,"_")[[1]][1:2],collapse =  "_")))
  
  subjectOutDT <- lapply(
    subjV,
    function(sampleID){
      #run both methods
      csSAMfitTemp <- csSAMcsFit(
        inputPeakDT[grep(sampleID,sample)],
        inputPropDT[grep(sampleID,sample)])
      
      #label them
      csSAMfitTemp[,method:='csSAM']
      
      return(
        csSAMfitTemp[,sampleID:=sampleID]
      )
    }
  )
  
  return(
    rbindlist(
      subjectOutDT
    )
  )
}


## ----------------------------------------------------------------------------------------------------------------------------------

adMixsingleSubjectDecompositionDT <- singleSampleDecRunner(mixPeakDT,mixPropDT)


## ----------------------------------------------------------------------------------------------------------------------------------
cmpDT <- admixture_input_DT[adMixsingleSubjectDecompositionDT,
                            o=.(status==cell,
                                subj==sampleID,
                                rn==peak)]


## ----------------------------------------------------------------------------------------------------------------------------------
library(ggpubr)
library(ggthemes)

cmpDT[,Age:='Old']
cmpDT[grep("Y",subj),Age:='Young']

cmpDT$Age <- factor(cmpDT$Age, levels = c("Young", "Old"))

cmpDT$status <- factor(cmpDT$status, levels = c("Activated", "Quiescent", "Neurons", "Oligodendrocytes"))

ggplot(cmpDT[method=='csSAM'],
       aes(x=value,
           y=mean
       )) + 
  geom_point(aes(color=status), alpha = 0.3) + 
  facet_grid( ~ status)+
  scale_color_manual(values = col.ls)+
  theme(legend.position = 'bottom') + 
  xlab('Measured lipid intensity') +
  ylab("Estimated lipid intensity")+
  theme_classic() +
  theme(axis.text = element_text(colour = "black", 
                                 face = "plain",
                                 size = 8))+
  stat_cor(label.sep = "\n", size = 2.5)+
  theme(legend.position = "none")
ggsave("./Figure_Panels/fig.S7A.pdf", width = 6, height = 2, useDingbats=FALSE)


## ----------------------------------------------------------------------------------------------------------------------------------
permute_prop_DT <- function(inputDT){
  
  nCells   <- dim(inputDT)[1]*(dim(inputDT)[2]-1)
  inputMat <- t(
    apply(
      matrix(
        sample(nCells),
        ncol=4
      ),1,
      function(i) i/sum(i)
    ))
  row.names(inputMat) <- inputDT$sample
  
  outDT  <- as.data.table(
    inputMat,
    keep.rownames = T)
  setnames(outDT,colnames(inputDT))
  
  return(outDT)
}


## ----------------------------------------------------------------------------------------------------------------------------------
set.seed(3333)

adMixPermutationsDT <- rbindlist(
  lapply(1:100,
         function(nP)
           singleSampleDecRunner(
             inputPeakDT = mixPeakDT,
             inputPropDT = permute_prop_DT(mixPropDT))[,nP:=nP]))


## ----------------------------------------------------------------------------------------------------------------------------------
permuationsCmpDT <- admixture_input_DT[adMixPermutationsDT,
                                       o=.(status==cell,
                                           subj==sampleID,
                                           rn==peak)]
permuationsCmpDT[,Age:='Old']
permuationsCmpDT[grep("Y",subj),Age:='Young']

permuationsCmpDT$Age <- factor(permuationsCmpDT$Age,
                               levels = c("Young", "Old"))


## ----------------------------------------------------------------------------------------------------------------------------------
permuationsCmpDT$status <- factor(
  permuationsCmpDT$status, levels = c("Activated", 
                                      "Quiescent", 
                                      "Neurons", 
                                      "Oligodendrocytes"))

cmpDT$status <- factor(cmpDT$status, levels = c("Activated", 
                                                "Quiescent", 
                                                "Neurons", 
                                                "Oligodendrocytes"))
ggplot(permuationsCmpDT[,cor(value,mean),
                        .(method,status,nP)][method=='csSAM'],
       aes(x=V1,y=..count.., fill=status)) + 
  
  geom_density(alpha = 0.75) + 
  
  scale_fill_manual(values = col.ls) + 
  
  xlab("Correlation coefficient between estimated and measured lipid intensity") +
  labs(fill='cell type')+
  
  theme(legend.position = 'none') +
  facet_grid(.~status) +
  
  geom_vline(data = cmpDT[method=='csSAM'][
    ,cor(value,mean),
    .(method,status)],
    aes(xintercept=V1,color=status),linetype=2) +
  scale_color_manual(values = col.ls) + 
  
  labs(fill='cell type')+
  labs(color='cell type') +
  
  theme_classic() + 
  theme(axis.text = element_text(colour = "black", 
                                 face = "plain",
                                 size = 8))+
  theme(legend.position = 'none')
ggsave("./Figure_Panels/fig.S7B.pdf", width = 6, height = 2, useDingbats=FALSE)

