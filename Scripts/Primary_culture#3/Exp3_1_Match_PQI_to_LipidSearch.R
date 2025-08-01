## Primary culture #3 data, matching Progenesis extracted peaks with LipidSearch annotation on peaks

# Match with the Lipid Molecule export from LipidSearch, use "TotalGrade" as a lipid quality metric

# Matching criteria:
#1. m/z within the rage of +/-0.01
#2. ppm within the rage of +/-10
#3. RT within the rage of +/-0.025*RT from LipidSearch annotated peaks

# Will also try matching with lipidID output, in addition to using lipidmolecule output
rm(list=ls())
setwd(rstudioapi::getActiveProject())
library(openxlsx)
library(tidyverse)
library('bitops')
## Process progensis data
### Positive mode
#==== Progenesis data organization - positive ======
PQI_pos <- read.csv("./Input_Data/231114_Xiaoai_Lipidomics_Batch3_pos_forR.csv", stringsAsFactors = F)
smp.pos <- PQI_pos %>% 
  select(., c("Compound",starts_with("XZ_"))) %>% 
  column_to_rownames(var = "Compound")
smp.pos[smp.pos == 0] <- NA
smp.pos <- smp.pos[-which(rowSums(is.na(smp.pos) == T) >= 0.5*ncol(smp.pos)),] #15798

trim.pos <- PQI_pos %>% 
  filter(Compound %in% rownames(smp.pos)) %>% 
  filter(Retention.time..min. > 0.5) 
smpl.trim <- trim.pos %>% 
  select(XZ_83:XZ_104) %>% 
  select(-c("XZ_86_PBS", "XZ_96_PBS"))
med_bio <- apply(smpl.trim,1,median, na.rm = TRUE)

trim.pos$med_bio <- med_bio

trim.pos <- trim.pos %>% 
  rowwise() %>% 
  mutate(., StN = med_bio/mean(XZ_86_PBS, XZ_96_PBS)) %>% 
  filter(StN >= 2.0) #6646

##sanity check
i = 5023
t.trim <- trim.pos %>% 
  select(XZ_83:XZ_104) %>% 
  select(-c("XZ_86_PBS", "XZ_96_PBS"))
trim.pos$med_bio[i] == apply(t.trim,1, median, na.rm = TRUE)[i]

#==== All Positive data from LipidSearch ======############
All.pos.mol <- read.csv("./Input_Data/Batch_3_LipidMolecule_POS_align.txt", sep = "\t")
# All.pos.ID <- read.csv("./Scripts/2023_dataset/Batch_3/LipidSearch_alignment/021124_Batch3_allsample_pos_align_LipidIon.txt", sep = "\t")

## add charge to get m/z
Charge <- c("22.989770",
            "1.007825",
            "18.034374",
            "-17.00274")

mz.tbl <- tibble(unique(All.pos.mol$MainIon), Charge) %>%
  rename("Main_ion" = `unique(All.pos.mol$MainIon)`) %>%
  rename("Z" = "Charge") %>%
  mutate_at("Z", as.numeric)
mz.tbl

All.pos.mz <- All.pos.mol %>%
  mutate_at("CalcMass", as.numeric) %>%
  rowwise() %>%
  mutate(CalcMz = CalcMass + mz.tbl$Z[mz.tbl$Main_ion == MainIon]) %>%
  relocate(c(LipidMolec, LipidMolecGroup, CalcMz, CalcMass, MainIon), .after = ID)

# Only keep necessary columns
Match_pos <- trim.pos %>% 
  select(., c("Compound", "m.z", "Retention.time..min.", starts_with("QC_"), starts_with("XZ_")))

Match_pos$LipidIon <- ""
Match_pos$Class <- ""
Match_pos$IonFormula <- ""
Match_pos$RT_LS <- ""
Match_pos$mz_LS <- ""
Match_pos$Grade_LS <- ""
Match_pos$Ion_LS <- ""
Match_pos$Med_LS_Int <- ""
#pos and neg dataframe is annotated matrix from LipidSearch on all samples
#Match_pos dataframe is Progenesis peaks to be matched to lipidsearch annotation

for (i in 1:dim(Match_pos)[1]) {
  # i=sample(1:5954, 1) #3118
  range <- All.pos.mz[which(All.pos.mz$CalcMz < (Match_pos$m.z[i]+0.01) & All.pos.mz$CalcMz > (Match_pos$m.z[i]-0.01)),] #range is subset from QC data

  my.cur.ppms <- rep(0,dim(range)[1])
  my.cur.RTs <- rep(0,dim(range)[1])
  for (j in 1:dim(range)[1]) {
    my.cur.ppms[j] <- 1e6*((Match_pos[i,]$m.z-range[j,]$CalcMz)/Match_pos[i,]$m.z)
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
  my.mz <- range$CalcMz[my.same]
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
  filter(!LipidIon == "") #2725

i = 4 #666 = 1 ->NP 666 = 3 ->P 857 = 3 -> NP 857 = 1 -> P
test.QC <- All.pos.mz %>% 
  mutate(., t.Name = paste0(LipidMolec,
                            substr(MainIon, 2, nchar(MainIon))))
Match_pos.s$Grade_LS[i] == test.QC[test.QC$t.Name == Match_pos.s$LipidIon[i], "TotalGrade"]

Match3_pos.mol <- Match_pos
save(Match3_pos.mol, file = "./Output_Data/Matched_Batch3_pos_by_LipidMol.Rdata")


#==== Progenesis data organization - negative ======
#===new code
rm(list = ls())
PQI_neg <- read.csv("./Input_Data/231114_Xiaoai_Lipidomics_Batch3_neg_forR.csv", stringsAsFactors = F)
smp.neg <- PQI_neg %>% 
  select(., c("Compound",starts_with("XZ_"))) %>% 
  column_to_rownames(var = "Compound")
smp.neg[smp.neg == 0] <- NA
smp.neg <- smp.neg[-which(rowSums(is.na(smp.neg) == T) >= 0.5*ncol(smp.neg)),] #12674

trim.neg <- PQI_neg %>% 
  filter(Compound %in% rownames(smp.neg)) %>% 
  filter(Retention.time..min. > 0.5) #12447

smpl.trim.neg <- trim.neg %>% 
  select(XZ_83:XZ_104) %>% 
  select(-c("XZ_86_PBS", "XZ_96_PBS"))

med_bio.n <- apply(smpl.trim.neg,1,median, na.rm = TRUE)
trim.neg$med_bio <- med_bio.n
trim.neg <- trim.neg %>% 
  rowwise() %>% 
  mutate(., StN = med_bio/mean(XZ_86_PBS, XZ_96_PBS)) %>% 
  filter(StN >= 2.0) #5679

i = 503
t.trim <- trim.neg %>% 
  select(XZ_83:XZ_104) %>% 
  select(-c("XZ_86_PBS", "XZ_96_PBS"))
trim.neg$med_bio[i] == apply(t.trim,1, median, na.rm = TRUE)[i]


#==== QC Negative data from LipidSearch ======
All.neg.mol <- read.csv("./Input_Data/Batch_3_LipidMolecule_NEG_align.txt", sep = "\t")

## add charge to get m/z
Charge <- c(
            "2.015650",
            "1.007825",
            "15.023475",
            "44.997655",
            "59.013304")

mz.tbl <- tibble(unique(All.neg.mol$MainIon), Charge) %>% 
  rename("Main_ion" = `unique(All.neg.mol$MainIon)`) %>% 
  rename("Z" = "Charge") %>% 
  mutate_at("Z", as.numeric)
mz.tbl

All.neg.mz <- All.neg.mol %>% 
  mutate_at("CalcMass", as.numeric) %>% 
  rowwise() %>%
  mutate(CalcMz = case_when(
    str_detect(MainIon, regex("\\+")) ~ CalcMass + mz.tbl$Z[mz.tbl$Main_ion == MainIon],
    str_detect(MainIon, regex("\\-")) ~ CalcMass - mz.tbl$Z[mz.tbl$Main_ion == MainIon]
  )) %>%
  relocate(c(LipidMolec, LipidMolecGroup, CalcMz, CalcMass, MainIon), .after = ID) 

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
  # i=sample(1:8000, 1) #3397
  range <- All.neg.mz[which(All.neg.mz$CalcMz < Match_neg$m.z[i]+0.01 & All.neg.mz$CalcMz > Match_neg$m.z[i]-0.01),] #range is subset from QC data
  my.cur.ppms <- rep(0,dim(range)[1])
  my.cur.RTs <- rep(0,dim(range)[1])
  for (j in 1:dim(range)[1]) {
    my.cur.ppms[j] <- 1e6*((Match_neg[i,]$m.z-range[j,]$CalcMz)/Match_neg[i,]$m.z)
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
  filter(!LipidIon == "") #1396

i = 5 #666 = 1 ->NP 666 = 3 ->P 857 = 3 -> NP 857 = 1 -> P
test.QC.n <- All.neg.mz %>% 
  mutate(., t.Name = paste0(LipidMolec,
                            substr(MainIon, 2, nchar(MainIon))))
Match_neg.s$Grade_LS[i] == test.QC.n[test.QC.n$t.Name == Match_neg.s$LipidIon[i], "TotalGrade"]

Match3_neg.mol <- Match_neg
save(Match3_neg.mol, file = "./Output_Data/Matched_Batch3_neg_by_LipidMol.Rdata")
