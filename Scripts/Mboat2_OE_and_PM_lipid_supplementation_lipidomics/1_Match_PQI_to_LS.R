## _Mboat2_ overexpression and young plasma membrane lipid supplementation lipidomics 
## Match Progenesis extracted peaks with LipidSearch annotation on peaks

# Steps:
#1. Keep peaks that are detected in at least half of the samples
#2. Remove peaks that are in the first 0.5min of the elution
#3. Peaks with Signal/Noise of >2

# Matching criteria:
#1. m/z within the rage of +/-0.01
#2. ppm within the rage of +/-10
#3. RT within the rage of +/-0.025*RT from LipidSearch annotated peaks

rm(list=ls())
setwd(rstudioapi::getActiveProject())
library(openxlsx)
library(tidyverse)
library('bitops')
## Process progensis data
### Positive mode
#==== Progenesis data organization - positive ======
PQI_pos <- read.csv("./Input_Data/230723_Xiaoai_Lipidomics_Batch2_pos_forR.csv", stringsAsFactors = F)
smp.pos <- PQI_pos %>% 
  select(., c("Compound",starts_with("XZ_"))) %>% 
  column_to_rownames(var = "Compound")
smp.pos[smp.pos == 0] <- NA
smp.pos <- smp.pos[-which(rowSums(is.na(smp.pos) == T) >= 0.5*ncol(smp.pos)),] #13926

trim.pos <- PQI_pos %>% 
  filter(Compound %in% rownames(smp.pos)) %>% 
  filter(Retention.time..min. > 0.5) #13536
smpl.trim <- trim.pos %>% 
  select(XZ_43:XZ_82)
med_bio <- apply(smpl.trim,1,median, na.rm = TRUE)

trim.pos$med_bio <- med_bio
trim.pos <- trim.pos %>% 
  rowwise() %>% 
  mutate(., StN = med_bio/mean(PBS_Prep_Blank_2, PBS_Prep_Blank_3)) %>% 
  filter(StN >= 2.0) #5954

##sanity check
i = 523
t.trim <- trim.pos %>% 
  select(XZ_43:XZ_82)
trim.pos$med_bio[i] == apply(t.trim,1, median, na.rm = TRUE)[i]

#==== QC Positive data from LipidSearch ======############
QC.pos.mol <- read.csv("./Input_Data/Batch_2_pos_5QC_9samples_LipidMolecule_export", sep = "\t")

## add charge to get m/z
Charge <- c("1.007825",
            "22.989770",
            "18.034374",
            "7.016005",
            "38.963708",
            "-17.00274",
            "46.065674 ")

mz.tbl <- tibble(unique(QC.pos.mol$MainIon), Charge) %>% 
  rename("Main_ion" = `unique(QC.pos.mol$MainIon)`) %>% 
  rename("Z" = "Charge") %>% 
  mutate_at("Z", as.numeric)
mz.tbl

QC.pos.mz <- QC.pos.mol %>% 
  mutate_at("CalcMass", as.numeric) %>% 
  rowwise() %>% 
  mutate(CalcMZ = CalcMass + mz.tbl$Z[mz.tbl$Main_ion == MainIon]) %>% 
  relocate(c(LipidMolec, LipidMolecGroup, CalcMZ, CalcMass, MainIon), .after = ID) 

## Check
t=789
d <- QC.pos.mz$CalcMZ[t] - QC.pos.mz$CalcMass[t]
all.equal(d, mz.tbl$Z[mz.tbl$Main_ion == QC.pos.mz$MainIon[t]])

# Only keep necessary columns
Match_pos <- trim.pos %>% 
  select(., c("Compound", "m.z", "Retention.time..min.", starts_with("QC_"), starts_with("XZ_"))) %>% 
  arrange(m.z)

Match_pos$LipidIon <- ""
Match_pos$Class <- ""
Match_pos$IonFormula <- ""
Match_pos$RT_LS <- ""
Match_pos$mz_LS <- ""
Match_pos$Grade_LS <- ""
Match_pos$Ion_LS <- ""
Match_pos$Med_LS_Int <- ""
#pos and neg dataframe is annotated matrix from LipidSearch on QC samples
#Match_pos dataframe is Progenesis peaks to be matched to QC

for (i in 1:dim(Match_pos)[1]) {
  # i=sample(1:5954, 1) #538 #3974
  #i = 288
  range <- QC.pos.mz[which(QC.pos.mz$CalcMZ < Match_pos$m.z[i]+0.01 & QC.pos.mz$CalcMZ > Match_pos$m.z[i]-0.01),] #range is subset from QC data
  my.cur.ppms <- rep(0,dim(range)[1])
  my.cur.RTs <- rep(0,dim(range)[1])
  for (j in 1:dim(range)[1]) {
    my.cur.ppms[j] <- 1e6*((Match_pos[i,]$m.z-range[j,]$CalcMZ)/Match_pos[i,]$m.z)
    my.cur.RTs[j] <- range$BaseRt[j]-Match_pos$Retention.time..min.[i]
  }
  my.same.mass <- which(bitAnd(my.cur.ppms < 10, my.cur.ppms > -10) > 0)
  my.same.RT <- which(bitAnd(my.cur.RTs < Match_pos$Retention.time..min.[i]*0.025, my.cur.RTs > -Match_pos$Retention.time..min.[i]*0.025) > 0)
  my.same <- intersect(my.same.mass, my.same.RT)
  my.ID <- paste0(range$LipidMolec[my.same],
                  substr(range$MainIon[my.same], 2, nchar(range$MainIon[my.same]))) # add ion to LipidIon 
  my.Formula <- range$MolFormula[my.same]
  my.Class <- range$ClassKey[my.same]
  my.RT <- range$BaseRt[my.same]
  my.mz <- range$CalcMass[my.same]
  my.grade <- range$TotalGrade[my.same]
  my.Ion <- range$MainIon[my.same]
  my.Int <- range$MedArea.s1.[my.same]
  if(length(my.same) >= 1){
    Match_pos$LipidIon[i] <- paste(my.ID, collapse = "| ")
    Match_pos$Class[i] <- paste(as.character(my.Class), collapse = "| ")
    Match_pos$IonFormula[i] <- paste(as.character(my.Formula), collapse = "| ")
    Match_pos$RT_LS[i] <- paste(my.RT, collapse = "| ")
    Match_pos$mz_LS[i] <- paste(my.mz, collapse = "| ")
    Match_pos$Grade_LS[i] <- paste(my.grade, collapse = "| ")
    Match_pos$Ion_LS[i] <- paste(my.Ion, collapse = "| ")
    Match_pos$Med_LS_Int[i] <- paste(my.Int, collapse = "| ")
  }else{
  }
}

Match_pos$Mode <- "Accucore_C18_pos"

##sanity check
Match_pos.s <- Match_pos %>%
  relocate(LipidIon, .after = Compound) %>%
  filter(!LipidIon == "") #2906 with +/- 5ppm, #3136 with +/- 10ppm

i = 112 #666 = 1 ->NP 666 = 3 ->P 857 = 3 -> NP 857 = 1 -> P
test.QC <- QC.pos.mz %>% 
  mutate(., t.Name = paste0(LipidMolec,
                           substr(MainIon, 2, nchar(MainIon))))
Match_pos.s$Grade_LS[i] == unlist(unique(test.QC[test.QC$t.Name == Match_pos.s$LipidIon[i], "TotalGrade"]))

Match2_pos.mol <- Match_pos
save(Match2_pos.mol, file = "./Output_Data/Matched_Batch2_pos_by_LipidMol.Rdata")


#==== Progenesis data organization - negative ======
#===new code
rm(list = ls())
PQI_neg <- read.csv("./Input_Data/230723_Xiaoai_Lipidomics_Batch2_neg_forR.csv", stringsAsFactors = F)
smp.neg <- PQI_neg %>% 
  select(., c("Compound",starts_with("XZ_"))) %>% 
  column_to_rownames(var = "Compound")
smp.neg[smp.neg == 0] <- NA
smp.neg <- smp.neg[-which(rowSums(is.na(smp.neg) == T) >= 0.5*ncol(smp.neg)),] #10463

trim.neg <- PQI_neg %>% 
  filter(Compound %in% rownames(smp.neg)) %>% 
  filter(Retention.time..min. > 0.5) #10081

smpl.trim.neg <- trim.neg %>% 
  select(XZ_43:XZ_82)

med_bio.n <- apply(smpl.trim.neg,1,median, na.rm = TRUE)
trim.neg$med_bio <- med_bio.n
trim.neg <- trim.neg %>% 
  rowwise() %>% 
  mutate(., StN = med_bio/mean(PBS_Prep_Blank_3)) %>% 
  filter(StN >= 2.0) #5085

i = 503
t.trim <- trim.neg %>% 
  select(XZ_43:XZ_82)
trim.neg$med_bio[i] == apply(t.trim,1, median, na.rm = TRUE)[i]


#==== QC Negative data from LipidSearch ======

QC.neg.mol <- read.csv("./Input_Data/Batch_2_neg_3QC_9samples_LipidMolecule_export", sep = "\t")

## add charge to get m/z
Charge <- c("1.007825",
            "2.015650",
            "44.997655",
            "15.023475",
            "59.013304")

mz.tbl <- tibble(unique(QC.neg.mol$MainIon), Charge) %>% 
  rename("Main_ion" = `unique(QC.neg.mol$MainIon)`) %>% 
  rename("Z" = "Charge") %>% 
  mutate_at("Z", as.numeric)
mz.tbl

QC.neg.mz <- QC.neg.mol %>% 
  mutate_at("CalcMass", as.numeric) %>% 
  rowwise() %>%
  mutate(CalcMZ = case_when(
    str_detect(MainIon, regex("\\+")) ~ CalcMass + mz.tbl$Z[mz.tbl$Main_ion == MainIon],
    str_detect(MainIon, regex("\\-")) ~ CalcMass - mz.tbl$Z[mz.tbl$Main_ion == MainIon]
  )) %>%
  relocate(c(LipidMolec, LipidMolecGroup, CalcMZ, CalcMass, MainIon), .after = ID) 

# Only keep necessary columns
Match_neg <- trim.neg %>% 
  select(., c("Compound", "m.z", "Retention.time..min.", starts_with("QC_"), starts_with("XZ_")))

Match_neg$LipidIon <- ""
Match_neg$Class <- ""
Match_neg$IonFormula <- ""
Match_neg$RT_LS <- ""
Match_neg$mz_LS <- ""
Match_neg$Grade_LS <- ""
Match_neg$Ion_LS <- ""
Match_neg$Med_LS_Int <- ""
#pos and neg dataframe is annotated matrix from LipidSearch on QC samples
#Match_pos dataframe is Progenesis peaks to be matched to QC

for (i in 1:dim(Match_neg)[1]) {
  # i=sample(1:18000, 1) #538
  range <- QC.neg.mz[which(QC.neg.mz$CalcMZ < Match_neg$m.z[i]+0.01 & QC.neg.mz$CalcMZ > Match_neg$m.z[i]-0.01),] #range is subset from QC data
  my.cur.ppms <- rep(0,dim(range)[1])
  my.cur.RTs <- rep(0,dim(range)[1])
  for (j in 1:dim(range)[1]) {
    my.cur.ppms[j] <- 1e6*((Match_neg[i,]$m.z-range[j,]$CalcMZ)/Match_neg[i,]$m.z)
    my.cur.RTs[j] <- range$BaseRt[j]-Match_neg$Retention.time..min.[i]
  }
  my.same.mass <- which(bitAnd(my.cur.ppms < 10, my.cur.ppms > -10) > 0)
  my.same.RT <- which(bitAnd(my.cur.RTs < Match_neg$Retention.time..min.[i]*0.025, my.cur.RTs > -Match_neg$Retention.time..min.[i]*0.025) > 0)
  my.same <- intersect(my.same.mass, my.same.RT)
  my.ID <- paste0(range$LipidMolec[my.same],
                  substr(range$MainIon[my.same], 2, nchar(range$MainIon[my.same]))) # add ion to LipidIon 
  my.Formula <- range$MolFormula[my.same]
  my.Class <- range$ClassKey[my.same]
  my.RT <- range$BaseRt[my.same]
  my.mz <- range$CalcMass[my.same]
  my.grade <- range$TotalGrade[my.same]
  my.Ion <- range$MainIon[my.same]
  my.Int <- range$MedArea.s1.[my.same]
  if(length(my.same) >= 1){
    Match_neg$LipidIon[i] <- paste(my.ID, collapse = "| ")
    Match_neg$Class[i] <- paste(as.character(my.Class), collapse = "| ")
    Match_neg$IonFormula[i] <- paste(as.character(my.Formula), collapse = "| ")
    Match_neg$RT_LS[i] <- paste(my.RT, collapse = "| ")
    Match_neg$mz_LS[i] <- paste(my.mz, collapse = "| ")
    Match_neg$Grade_LS[i] <- paste(my.grade, collapse = "| ")
    Match_neg$Ion_LS[i] <- paste(my.Ion, collapse = "| ")
    Match_neg$Med_LS_Int[i] <- paste(my.Int, collapse = "| ")
  }else{
  }
}

Match_neg$Mode <- "Accucore_C18_neg"

##sanity check
Match_neg.s <- Match_neg %>%
  relocate(LipidIon, .after = Compound) %>%
  filter(!LipidIon == "") #1673

i = 8 #666 = 1 ->NP 666 = 3 ->P 857 = 3 -> NP 857 = 1 -> P
test.QC.n <- QC.neg.mz %>% 
  mutate(., t.Name = paste0(LipidMolec,
                            substr(MainIon, 2, nchar(MainIon))))
Match_neg.s$Grade_LS[i] == test.QC.n[test.QC.n$t.Name == Match_neg.s$LipidIon[i], "TotalGrade"]

Match2_neg.mol <- Match_neg
save(Match2_neg.mol, file = "./Output_Data/Matched_Batch2_neg_by_LipidMol.Rdata")
