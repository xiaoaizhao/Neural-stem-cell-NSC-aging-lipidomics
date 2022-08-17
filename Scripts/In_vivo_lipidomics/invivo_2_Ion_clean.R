##Second processing script for in vivo sorted NSC lipidomics

##Steps:

##1. Ion clean for In vivo sorted qNSC data
##2. Manually remove one ion as it was identified as a spike-in in previous script
##3. Manually add cholesterol annotation

#Result:
#159 lipids total with correct ion.
rm(list=ls())
library(tidyverse)
setwd(rstudioapi::getActiveProject())

load(file = "./Output_Data/Multi_ID_corrected_invivo.Rdata")

##filter to get annotated lipids
df.id <- exvivo.df %>%
  filter(., !LipidIon == "0") %>% #313
  filter(., !LipidIon.FA == "Cer(m32:0)+NH4") #312, in script (invivo_0_Label_IS.R), this ion was identified as Ceramide (d7) IS, so remove from matrix

##Manually label endogenous Cholesterol, labelled as "ChE(0:0)+H-H2O", with RT at around 10mins
#Chol
Chol <- exvivo.df %>%
  filter(., m.z < 369.4 & m.z >369.3 & Mode == "pos"& Retention.time..min.>9 & Retention.time..min.<10) %>%
  mutate(., LipidIon.FA = "Cholesterol(0:0)+H-H2O") %>%
  mutate(., Class = "Cholesterol")

All.id <- bind_rows(df.id, Chol) 


##import ion list to keep for each specific lipid, this was determined from previous experiment based on the most abundant ion of each lipid class
ion.tbl <- read.csv("./Input_Data/Final_ion_list_for_cleanUp_052620_for2017LC.csv", stringsAsFactors = F)

##add ion info to each lipid
Ion.clean.159 <-All.id %>%
  mutate(., Ion = paste0("[", substr(LipidIon.FA, str_locate(LipidIon.FA, "\\)")+1, nchar(LipidIon.FA)), "]")) %>%
  group_by(Compound) %>%
  group_modify(~ {
    .x %>%
      mutate(., ion_match = ifelse(Ion == ion.tbl$Ion[match(Class, ion.tbl$Class)], "T", "F"))
  }) %>%
  filter(., ion_match == "T")  %>% 
  arrange(., Ion) ##159 after ion cleanup.

save(Ion.clean.159, file = paste0("./Output_Data/Invivo_IonClean_final_159_lipids.Rdata"))

##Inspect what ion was used to subset each lipid class
unique(Ion.clean.159$Class[Ion.clean.159$Ion == "[-H]"])
# [1] "LPE" "LPI" "PI"  "PS" 
unique(Ion.clean.159$Class[Ion.clean.159$Ion == "[+H]"])
# [1] "AcCa" "LPC"  "PC"   "PE"   "SM"
unique(Ion.clean.159$Class[Ion.clean.159$Ion == "[+NH4]"])
# [1] "TG"
unique(Ion.clean.159$Class[Ion.clean.159$Ion == "[+H-H2O]"])
# [1] "Cer"         "Cholesterol" "Hex1Cer"
