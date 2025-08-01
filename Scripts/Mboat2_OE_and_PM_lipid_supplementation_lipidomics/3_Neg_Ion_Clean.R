## Ion clean for each class
## Remove multiple annotation from negative mode data and clean based on ions used for quantification for each class.

## Steps (1673 ion with annotation)
#1. lipid with unique ID (734), 226 left after ion clean.
#2. lipid with multiple ID (939) - general rules: keep lipid that has the highest overall grade
  #1. lipids with one identification that has the highest grade (257, m.mult.1)
  #2. lipids with multiple identification with identical highest grade (682), then take lipids with the highest median intensity
    #1. 668 unique lipid can be identified based on the highest median intensity, m.mlt.grade
    #2. 14 lipids have identical value of highest median intensity, randomly pick one identification(first one in the list), m.mlt.grade.2c
# final cleaned positive ion = 496
rm(list=ls())
setwd(rstudioapi::getActiveProject())
library(openxlsx)
library(tidyverse)

load("./Output_Data/Matched_Batch2_neg_by_LipidMol.Rdata") #5058
Match_neg.s <- Match2_neg.mol %>%
  relocate(LipidIon, .after = Compound) %>%
  filter(!LipidIon == "") #1673

ion.tbl <- read.csv("./Input_Data/Final_ion_list_for_cleanUp_052620_for2017LC.csv", stringsAsFactors = F)

#==== ion clean on lipids with one unique identification ====

uni.ID <- Match_neg.s %>% 
  filter(., !grepl("\\|", LipidIon)) #734 lipids with only one annotation

uni.ID.c <-uni.ID %>%
  mutate(., Class1 = Class) %>% 
  mutate(., Ion = paste0("[", substr(LipidIon, str_locate(LipidIon, "\\)")+1, nchar(LipidIon)), "]")) %>%
  group_by(Class1) %>%
  group_modify(~ {
    .x %>%
      mutate(., ion_match = ifelse(Ion == ion.tbl$Ion[match(Class, ion.tbl$Class)], "T", "F"))
  }) %>%
  relocate(c("Ion","ion_match", "Grade_LS"), .after =  "LipidIon") %>% 
  filter(., ion_match == "T")  %>%
  arrange(., Ion) ##226 after ion cleanup.

uni.ID.c$RT_LS <- as.numeric(uni.ID.c$RT_LS)
uni.ID.c$mz_LS <- as.numeric(uni.ID.c$mz_LS)
uni.ID.c$Med_LS_Int <- as.numeric(uni.ID.c$Med_LS_Int)

#==== correct multiple ID by detection grade ====
## pick features that has the highest overall grading
multi.id <- Match_neg.s %>% 
  filter(., grepl("\\|", LipidIon)) %>% #939
  mutate(., IDs = str_split(LipidIon, "\\| ")) %>%
  mutate(., Grade = str_split(Grade_LS, "\\| ")) %>%
  mutate(., RT_all = str_split(RT_LS, "\\| ")) %>%
  mutate(., mz_all = str_split(mz_LS, "\\| ")) %>%
  mutate(., Int_all = str_split(Med_LS_Int, "\\| ")) %>% 
  relocate(c("IDs", "Grade", "RT_all"),
           .after =  "LipidIon") %>%
  rowwise() %>% 
  mutate(., ID.toKeep = list(IDs[which(unlist(Grade) == min(unlist(Grade)))])) %>%
  mutate(., RT.toKeep = list(RT_all[which(unlist(Grade) == min(unlist(Grade)))])) %>% 
  mutate(., MZ.toKeep = list(mz_all[which(unlist(Grade) == min(unlist(Grade)))])) %>% 
  mutate(., MedInt.toKeep = list(Int_all[which(unlist(Grade) == min(unlist(Grade)))])) %>% 
  mutate(., Grade.toKeep = list(Grade[which(unlist(Grade) == min(unlist(Grade)))])) %>% 
  arrange(desc(Compound)) %>% 
  relocate(c("ID.toKeep","RT.toKeep", "IDs", "Grade"),
           .after =  "LipidIon")

# these lipid have one single identification that has the highest confidence
## DF1 of lipids that are matched with multiple ID, pick the one ID with the highest identification grade 
m.mult.1 <- multi.id %>% 
  filter(length(unlist(ID.toKeep)) == 1) #257
m.mult.1$ID.toKeep <- as.character(m.mult.1$ID.toKeep)
m.mult.1$MZ.toKeep <- as.numeric(m.mult.1$MZ.toKeep)
m.mult.1$RT.toKeep <- as.numeric(m.mult.1$RT.toKeep)
m.mult.1$MedInt.toKeep <- as.numeric(m.mult.1$MedInt.toKeep)
m.mult.1$Grade.toKeep <- as.character(m.mult.1$Grade.toKeep)

# these lipid have multiple identifications that have similar level of confidence
m.mult.2 <- multi.id %>% 
  filter(!length(unlist(ID.toKeep)) == 1) #682

##theselipids have multiple highest grading, then do second round of ranking base on median intensity from LS
## all return unique identification with this approach
m.mlt.grade <- m.mult.2 %>% #682
  rowwise() %>% 
  mutate(., ID.toKeepI = list(ID.toKeep[which(unlist(MedInt.toKeep) == max(as.numeric(unlist(MedInt.toKeep))))])) %>%
  mutate(., RT.toKeepI = list(RT.toKeep[which(unlist(MedInt.toKeep) == max(as.numeric(unlist(MedInt.toKeep))))])) %>% 
  mutate(., MZ.toKeepI = list(MZ.toKeep[which(unlist(MedInt.toKeep) == max(as.numeric(unlist(MedInt.toKeep))))])) %>% 
  mutate(., Grade.toKeepI = list(Grade.toKeep[which(unlist(MedInt.toKeep) == max(as.numeric(unlist(MedInt.toKeep))))])) %>% 
  mutate(., MedInt.toKeepI = list(MedInt.toKeep[which(unlist(MedInt.toKeep) == max(as.numeric(unlist(MedInt.toKeep))))])) %>% 
  relocate(c("ID.toKeepI","MedInt.toKeepI", "MedInt.toKeep"),
           .after =  "ID.toKeep")  %>% 
  filter(length(unlist(ID.toKeepI)) == 1) #668

m.mlt.grade2 <-  m.mult.2 %>% 
  rowwise() %>% 
  mutate(., ID.toKeepI = list(ID.toKeep[which(unlist(MedInt.toKeep) == max(as.numeric(unlist(MedInt.toKeep))))])) %>%
  mutate(., RT.toKeepI = list(RT.toKeep[which(unlist(MedInt.toKeep) == max(as.numeric(unlist(MedInt.toKeep))))])) %>% 
  mutate(., MZ.toKeepI = list(MZ.toKeep[which(unlist(MedInt.toKeep) == max(as.numeric(unlist(MedInt.toKeep))))])) %>% 
  mutate(., Grade.toKeepI = list(Grade.toKeep[which(unlist(MedInt.toKeep) == max(as.numeric(unlist(MedInt.toKeep))))])) %>% 
  mutate(., MedInt.toKeepI = list(MedInt.toKeep[which(unlist(MedInt.toKeep) == max(as.numeric(unlist(MedInt.toKeep))))])) %>% 
  relocate(c("ID.toKeepI","MedInt.toKeepI", "MedInt.toKeep"),
           .after =  "ID.toKeep")  %>% 
  filter(!length(unlist(ID.toKeepI)) == 1) #14

m.mlt.grade.2c <- m.mlt.grade2 %>% 
  rowwise() %>% 
  mutate(., ID.toKeepII = unlist(ID.toKeepI)[1]) %>%
  mutate(., RT.toKeepII = unlist(RT.toKeepI)[1]) %>% 
  mutate(., MZ.toKeepII = unlist(MZ.toKeepI)[1]) %>% 
  mutate(., Grade.toKeepII = unlist(Grade.toKeepI)[1]) %>% 
  mutate(., MedInt.toKeepII = unlist(MedInt.toKeepI)[1]) %>% 
  relocate(c("ID.toKeepII","MedInt.toKeepII", "MedInt.toKeepI"),
           .after =  "ID.toKeep") 

rmv.I <- function(x){
  x <- ifelse(
    grepl("I$", x),
    str_replace_all(x, "I", ""),
    x
  )
}

m.df.1 <- m.mult.1
m.df.1$ID.toKeep <- as.character(m.df.1$ID.toKeep)
m.df.1$Grade.toKeep <- as.character(m.df.1$Grade.toKeep)
m.df.1$MedInt.toKeep <- as.numeric(m.df.1$MedInt.toKeep)
m.df.1$RT.toKeep <- as.numeric(m.df.1$RT.toKeep)
m.df.1$MZ.toKeep <- as.numeric(m.df.1$MZ.toKeep)


m.df.2 <- m.mlt.grade %>% 
  select(-c(Grade.toKeep, RT.toKeep, MZ.toKeep, ID.toKeep, MedInt.toKeep)) %>% 
  rename_at(vars(contains("toKeep", ignore.case = FALSE)), ~rmv.I(.)) %>% 
  rename("ID.toKeep" = "D.toKeep") %>% 
  rename("MedInt.toKeep" = "Mednt.toKeep") 

m.df.2$ID.toKeep <- as.character(m.df.2$ID.toKeep)
m.df.2$Grade.toKeep <- as.character(m.df.2$Grade.toKeep)
m.df.2$MedInt.toKeep <- as.numeric(m.df.2$MedInt.toKeep)
m.df.2$RT.toKeep <- as.numeric(m.df.2$RT.toKeep)
m.df.2$MZ.toKeep <- as.numeric(m.df.2$MZ.toKeep)


m.df.3 <- m.mlt.grade.2c %>% 
  select(-c(Grade.toKeep, Grade.toKeepI,
            MedInt.toKeep, MedInt.toKeepI,
            RT.toKeep, RT.toKeepI,
            MZ.toKeep, MZ.toKeepI,
            ID.toKeep, ID.toKeepI)) %>% 
  rename_at(vars(contains("toKeep", ignore.case = FALSE)), ~rmv.I(.)) %>% 
  rename("ID.toKeep" = "D.toKeep") %>% 
  rename("MedInt.toKeep" = "Mednt.toKeep") 

m.df.3$ID.toKeep <- as.character(m.df.3$ID.toKeep)
m.df.3$Grade.toKeep <- as.character(m.df.3$Grade.toKeep)
m.df.3$MedInt.toKeep <- as.numeric(m.df.3$MedInt.toKeep)
m.df.3$RT.toKeep <- as.numeric(m.df.3$RT.toKeep)
m.df.3$MZ.toKeep <- as.numeric(m.df.3$MZ.toKeep)

Mult.to.clean <- bind_rows(m.df.1, m.df.2, m.df.3) %>% 
  select(-c(RT_LS, mz_LS, Grade_LS, Med_LS_Int, IDs, Grade, Int_all, mz_all, RT_all)) %>% 
  mutate(., LipidIon = ID.toKeep) %>% 
  mutate(., RT_LS = RT.toKeep) %>% 
  mutate(., mz_LS = MZ.toKeep) %>% 
  mutate(., Med_LS_Int = MedInt.toKeep) %>% 
  mutate(., Grade_LS = Grade.toKeep) %>% 
  select(-contains("toKeep")) %>% 
  rowwise() %>% 
  mutate(., Class = ifelse(LipidIon == "Ch+H-H2O",
                           "Cholesterol",
                           substr(LipidIon, 1, str_locate(LipidIon, "\\(")-1))) %>% 
  mutate(., Ion = ifelse(LipidIon == "Ch+H-H2O",
                         "[+H-H2O]",
                         paste0("[", substr(LipidIon, str_locate(LipidIon, "\\)")+1, nchar(LipidIon)), "]"))) 

Mlt.clean <- Mult.to.clean %>% #939
  mutate(., Class1 = Class) %>% 
  group_by(Class1) %>%
  group_modify(~ {
    .x %>%
      mutate(., ion_match = ifelse(Ion == ion.tbl$Ion[match(Class, ion.tbl$Class)], "T", "F"))
  }) %>%
  relocate(c("Ion","ion_match", "Grade_LS"), .after =  "LipidIon") %>% 
  filter(., ion_match == "T")  %>%
  arrange(., Ion) ##270 after ion cleanup.

M2PM.neg.clean <- bind_rows(uni.ID.c, Mlt.clean) #496

save(M2PM.neg.clean, file = "./Output_Data/M2PM_neg_ion_clean.Rdata")


### ===== manually inspect to see if the script is working as intended, i.e. filter based on highest grade, then by highest median intensity
## all checks out!!
i = 309
ins.df <- Match_neg.s[i, c("Compound", "LipidIon", "Grade_LS", "Med_LS_Int")]
ins.df

check.df <- M2PM.neg.clean[M2PM.neg.clean$Compound == ins.df$Compound,1:4]
check.df
