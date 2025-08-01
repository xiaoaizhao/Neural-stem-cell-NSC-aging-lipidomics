
###Filter to only keep lipids that are detected in at least half of the samples####=========================================================
##Convert 0 value to NA
##Count number of NAs, filter to only keep the number of NA < 1/2 of sample size
filter.func <- function(norm.df, smpl.size){
  Norm_filter <- norm.df
  Norm_filter[Norm_filter == 0] <- NA 
  Norm_filter <- Norm_filter %>%
    filter(., rowSums(is.na(Norm_filter)) < smpl.size *0.5)
}

###Filter to only keep lipids that are detected in at least half of the samples####
##same filter function, this one specifically for lipidyzer data, since the data matrix already has missing value presented as NA
filter.Ldz <- function(norm.df, smpl.size){
  Norm_filter <- norm.df %>%
    filter(., rowSums(is.na(norm.df)) < smpl.size *0.5)
  return(Norm_filter)
}

####Imputation function to replace 0 value in data matrix####================================================================================
##Input matrix has to be log2 transformed
##Convert all 0 value to NA, before log2 transformation
##Loop through each row (lipids), replace all NA with random value selected from 10% value of given lipids
impt.func <- function(Mednorm.df){
  Norm_0_replace <- Mednorm.df
  Norm_0_replace[Norm_0_replace == 0] <- NA ##convert all 0 value to 1 in order to do log2 transformation
  log2df <- log2(Norm_0_replace)
  imptdf <- log2df
  
  set.seed(12345)
  for (i in 1:nrow(imptdf)){
    b <- imptdf[i, which(is.na(imptdf[i,])==T)]
    imptdf[i, which(is.na(imptdf[i,])==T)]<- rnorm(n = length(b),
                                                   mean = quantile(as.numeric(imptdf[i,]),na.rm = T, 0.05),
                                                   sd = abs(0.3*quantile(as.numeric(imptdf[i,]),na.rm = T, 0.05))) 
  }
  return(imptdf)
}

####Imputation function to replace 0 value in data matrix####
##same imputation function, this one specifically for lipidyzer and Primary culture #2 data (Exp2), since the data matrix already has missing value presented as NA
impt.Ldz <- function(Mednorm.df){
  imptdf <- log2(Mednorm.df)
  set.seed(12345)
  for (i in 1:nrow(imptdf)){
    b <- imptdf[i, which(is.na(imptdf[i,])==T)]
    imptdf[i, which(is.na(imptdf[i,])==T)]<- rnorm(n = length(b),
                                                   mean = quantile(as.numeric(imptdf[i,]),na.rm = T, 0.05),
                                                   sd = abs(0.3*quantile(as.numeric(imptdf[i,]),na.rm = T, 0.05))) 
  }
  return(imptdf)
}



####Check the number of values that were replace by imputation = number of missing value to begin with####========================================
##Get number of missing value
##Get number of replaced value before and after imputation
##Make sure the number match
impt.check <- function(Mednorm.df, impt.df){
  Norm_0_replace <- Mednorm.df
  Norm_0_replace[Norm_0_replace == 0] = 1 ##convert all 0 value to 1 in order to do log2 transformation
  log2df <- log2(Norm_0_replace)
  n.missval = length(which(log2df==0)) 
  
  b <- setdiff(impt.df, log2df)
  c <- setdiff(log2df, impt.df)
  return(length(which(!b==c)) == n.missval)
}


####Check the number of values that were replace by imputation = number of missing value to begin with
##same function as above, this one specifically for lipidyzer and Primary culture #2 data (Exp2), since the data matrix already has missing value presented as NA
impt.check.Ldz <- function(Mednorm.df, impt.df){
  n.missval = length(Mednorm.df[is.na(Mednorm.df)])
  log2df <- log2(Mednorm.df)
  b <- setdiff(impt.df, log2df)
  c <- setdiff(log2df, impt.df)
  d <- as_tibble(!b==c)
  return(length(d[is.na(d)]) == n.missval)
}

####Check the number of values that were replace by imputation = number of missing value to begin with
##same function as above, this one specifically for z score data matrix used in Pct Mol analysis across data set. Since the data is alread scaled (by z score), there is no need for log2 transformation.
impt.check.zScore.mtx <- function(Mednorm.df, impt.df){
  n.missval = length(Mednorm.df[is.na(Mednorm.df)])
  b <- setdiff(impt.df, Mednorm.df)
  c <- setdiff(Mednorm.df, impt.df)
  d <- as_tibble(!b==c)
  return(length(d[is.na(d)]) == n.missval)
}


####Funciton for aggregating lipid intensity with any given double bond in each class=========================================================
##1. Subset by presence or absence of each given double bond, regardless of how many times (side chains) it was identified
##2. Tally total intensity with every given double bond number in each class.
db.tally <- function(df, value, sample){
Db <- paste0(":", rep(0:6))
  df_by_DB = list()
  sum_by_DB = list()
  for (numDB in Db) {
    df_by_DB[[numDB]] <- {{df}} %>%
      filter(., str_detect(SideChain, numDB)) %>%
      mutate(., DB_num = numDB) %>%
      group_by({{sample}}, Class) 
    sum_by_DB[[numDB]]  <- df_by_DB[[numDB]] %>%
      summarise(., Sum_DB = sum({{value}})) %>%
      mutate(., DB_num = numDB)
  }
  return(sum_by_DB)
}

db.tally.list <- function(df, value, sample){
  Db <- paste0(":", rep(0:6))
  df_by_DB = list()
  sum_by_DB = list()
  for (numDB in Db) {
    df_by_DB[[numDB]] <- {{df}} %>%
      filter(., any(str_detect(unlist(SideChain), numDB))) %>%
      mutate(., DB_num = numDB) %>%
      group_by({{sample}}, Class) 
    sum_by_DB[[numDB]]  <- df_by_DB[[numDB]] %>%
      summarise(., Sum_DB = sum({{value}})) %>%
      mutate(., DB_num = numDB)
  }
  return(sum_by_DB)
}
####Fold change on doubel bond composition between age####
DB.FC <- function(df){
  
  OY_mean <- {{df}} %>%
    group_by(Class, DB_num, Age) %>%
    summarise(., DB_Age_avg = mean(DB_Pct)) 
  
  OY_FC <- OY_mean %>%
    group_by(Class, DB_num) %>%
    group_modify(~ {
      .x %>%
        mutate(., log2FC_OvY = log2(DB_Age_avg/DB_Age_avg[Age == "Young"]))
    }) %>%
    filter(., log2FC_OvY !=0) %>%
    select(., -matches("Age|Age_Avg"))
  return(OY_FC)
}

DB.sat.FC <- function(df){
  
  OY_mean <- {{df}} %>%
    group_by(Class, Sat, Age) %>%
    summarise(., DB_Age_avg = mean(DB_Sat_Pct)) 
  
  OY_FC <- OY_mean %>%
    group_by(Class, Sat) %>%
    group_modify(~ {
      .x %>%
        mutate(., log2FC_OvY = log2(DB_Age_avg/DB_Age_avg[Age == "Young"]))
    }) %>%
    filter(., log2FC_OvY !=0) %>%
    select(., -matches("Age|Age_Avg"))
  return(OY_FC)
}


lpd.by.db <- function(df, value, sample){
  Db <- paste0(":", rep(0:6))
  df_by_DB = list()
  sum_by_DB = list()
  for (numDB in Db) {
    df_by_DB[[numDB]] <- {{df}} %>%
      filter(., str_detect(SideChain, numDB)) %>%
      mutate(., DB_num = numDB) %>%
      group_by({{sample}}, Class) 
    sum_by_DB[[numDB]]  <- df_by_DB[[numDB]] %>%
      mutate(., DB_num = numDB)
  }
  return(sum_by_DB)
}