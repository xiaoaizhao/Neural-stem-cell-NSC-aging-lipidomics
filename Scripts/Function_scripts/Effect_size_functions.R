#compute effect size [es] based on Hedge's g effect size.
#old vs young comparison 
#(pos ES -> higher in old)
#(neg ES -> higher in young)
#adapted from MetaIntegrator (CRAN)

## Effect size between old and young ####
es.g.func <- function(df, variable, Grp, value, smpl){
  smpl.size <- df %>%
    group_by({{Grp}}) %>%
    mutate(., n = length(unique({{smpl}}))) %>%
    summarise(., n.size = mean(n))

  mean.sd.df <- df %>%
    group_by({{variable}}, {{Grp}}) %>%
    summarize(Age_mean = mean({{value}}), Age_SD = sd({{value}})) %>%
    ungroup()
  
  eff.size.df <- left_join(mean.sd.df, smpl.size) %>%
    pivot_wider(names_from = {{Grp}}, values_from = c("Age_mean", "Age_SD", "n.size")) %>%
    mutate(., diff = Age_mean_Old - Age_mean_Young) %>%
    mutate(., sp = sqrt( ( (n.size_Old-1)*Age_SD_Old^2 + (n.size_Young-1)*Age_SD_Young^2 )/( n.size_Old + n.size_Young - 2 ) )) %>% #pooled SD
    mutate(., cf = 1 - 3/( 4*(n.size_Old + n.size_Young) - 9 )) %>% #correction factor for small sample size (n<50)
    mutate(., es_g = cf * diff/sp) %>% #Hedge's g effect size
    mutate(., se_g = sqrt( (n.size_Old+n.size_Young)/(n.size_Old*n.size_Young) + 0.5*es_g^2 /(n.size_Old+n.size_Young-3.94))) #Standard Error
}

## Effect size between quiescent and activated cells####
#(pos ES -> higher in Activated cell)
#(neg ES -> higher in Quiescent cell)
es.g.func.celltype <- function(df, variable, Grp, value, smpl){
  smpl.size <- df %>%
    group_by({{Grp}}) %>%
    mutate(., n = length(unique({{smpl}}))) %>%
    summarise(., n.size = mean(n))
  
  mean.sd.df <- df %>%
    group_by({{variable}}, {{Grp}}) %>%
    summarize(Cell_mean = mean({{value}}), Cell_SD = sd({{value}})) %>%
    ungroup()
  
  eff.size.df <- left_join(mean.sd.df, smpl.size) %>%
    pivot_wider(names_from = {{Grp}}, values_from = c("Cell_mean", "Cell_SD", "n.size")) %>%
    mutate(., diff = Cell_mean_Activated - Cell_mean_Quiescent ) %>%
    mutate(., sp = sqrt( ( (n.size_Quiescent-1)*Cell_SD_Quiescent^2 + (n.size_Activated-1)*Cell_SD_Activated^2 )/( n.size_Quiescent + n.size_Activated - 2 ) )) %>%
    mutate(., cf = 1 - 3/( 4*(n.size_Quiescent + n.size_Activated) - 9 )) %>%
    mutate(., es_g = cf * diff/sp) %>%
    mutate(., se_g = sqrt( (n.size_Quiescent+n.size_Activated)/(n.size_Quiescent*n.size_Activated) + 0.5*es_g^2 /(n.size_Quiescent+n.size_Activated-3.94)))
}

## Effect size between two different KOs####
#(pos ES -> higher in KO1)
#(neg ES -> higher in KO2)
es.g.func.KO<- function(df, variable, Grp, value, smpl){
  smpl.size <- df %>%
    group_by({{Grp}}) %>%
    mutate(., n = length(unique({{smpl}}))) %>%
    summarise(., n.size = mean(n))
  
  mean.sd.df <- df %>%
    group_by({{variable}}, {{Grp}}) %>%
    summarize(KO_mean = mean({{value}}), KO_SD = sd({{value}})) %>%
    ungroup()
  
  KO1 <- levels(df$KO)[1]
  KO2 <- levels(df$KO)[2]
  
  KO1_mean <- paste0("KO_mean_", KO1)
  KO2_mean <- paste0("KO_mean_", KO2)
  
  KO1_SD <- paste0("KO_SD_", KO1)
  KO2_SD <- paste0("KO_SD_", KO2)
  
  KO1.n.size <- paste0("n.size_", KO1)
  KO2.n.size <- paste0("n.size_", KO2)
  
  eff.size.df <- left_join(mean.sd.df, smpl.size) %>%
    pivot_wider(names_from = {{Grp}}, values_from = c("KO_mean", "KO_SD", "n.size")) %>%
    mutate(., diff = .data[[KO1_mean]] - .data[[KO2_mean]]) %>% 
    mutate(., sp = sqrt( ( (.data[[KO1.n.size]]-1)*.data[[KO1_SD]]^2 + (.data[[KO2.n.size]]-1)*.data[[KO2_SD]]^2 )/(.data[[KO1.n.size]] + .data[[KO2.n.size]] - 2 ) )) %>%
    mutate(., cf = 1 - 3/( 4*(.data[[KO1.n.size]] + .data[[KO2.n.size]]) - 9 )) %>%
    mutate(., es_g = cf * diff/sp) %>%
    mutate(., se_g = sqrt( (.data[[KO1.n.size]]+.data[[KO2.n.size]])/(.data[[KO1.n.size]]*.data[[KO2.n.size]]) + 0.5*es_g^2 /(.data[[KO1.n.size]]+.data[[KO2.n.size]]-3.94)))
}

###Up with age####
#Filter variables that has effect size >0 (up with age), create a ranked list based on effect size (absolute value high to low)
up.ES.DB.func <- function(df, met1, val1 ){
  up.org <- {{df}} %>%
    select(., c({{met1}}, {{val1}})) %>%
    arrange(., desc({{val1}})) %>%
    rowid_to_column(., var = "Rank_pos") %>%
    filter(., {{val1}} > 0) 
  up.org <- up.org %>%
    mutate(., Rank_Pct = Rank_pos/nrow(up.org)*100) %>%
    rename(., Effect_Size = es_g)
}

###Down with age####
#Filter variables that has effect size <0 (down with age), create a ranked list based on effect size (absolute value high to low)
down.ES.DB.func <- function(df, met1, val1 ){
  dw.org <- {{df}} %>%
    select(., c({{met1}}, {{val1}})) %>%
    arrange(., {{val1}}) %>%
    rowid_to_column(., var = "Rank_pos") %>%
    filter(., {{val1}} < 0) 
  dw.org <- dw.org %>%
    mutate(., Rank_Pct = Rank_pos/nrow(dw.org)*100) %>%
    rename(., Effect_Size = es_g)
}

####Append m.z. to effect size table####
add.mz <- function(dfw.mz, mzcolname, lipidcolname, es.df){
  noslash <- function(x, na.rm = FALSE){str_replace_all(x, "/", "_")}
  df <- {{dfw.mz}} %>%
    ungroup() %>%
    rename(., unfixedLipidIon = {{lipidcolname}}) %>%
    select(., c({{mzcolname}}, unfixedLipidIon)) %>%
    mutate(across(unfixedLipidIon),  noslash(unfixedLipidIon)) %>%
    rename(., LipidIonFix = `noslash(unfixedLipidIon)`)
  es.df <- es.df %>%
    ungroup() %>%
    rename(., LipidIonFix = {{lipidcolname}})
  es.w.mz <- left_join(es.df, df, by = "LipidIonFix")
}

####Convert all ion to -H to overlap with DESI data####
Conv.NegH <- function(df){
  ion <- c("-H", "+H", "+H-H2O", "+NH4")
  delta.mass <- c(1.00782504, -1.00782504, -1.00782504+18.010565, -18.034374) 
  ion.mass <- tibble(ion,delta.mass)

  all.mz <- df%>%
    mutate(., ion = substr(LipidIon, str_locate(LipidIon, "\\)")+1, nchar(LipidIon))) 
  all.to_negH <- left_join(all.mz, ion.mass, by = "ion") %>%
    mutate(., Original.mz = MZ_avg + delta.mass) %>%
    mutate(., DESI_m.z = Original.mz - 1.00782504) %>%
    mutate(., `del_-H` = MZ_avg- DESI_m.z) %>%
    mutate(., peak = round(DESI_m.z, 1))
}

####Reformating lipid identification string to overlap with dataset####
Lbr <- function(x, na.rm = FALSE){str_replace_all(x, "\\(", "\\\\(")}
Rbr <- function(x, na.rm = FALSE){str_replace_all(x, "\\)", "\\\\)")}


####Wilcox stat test and fold change calculation####
wilcox_stat <- function(df, var, grp){
  OY_stat <- {{df}} %>%
    group_by({{grp}}) %>%
    summarize(Wilcox_Pval = wilcox.test({{var}} ~ Age)$p.value) %>%
    mutate(., Padj = p.adjust(Wilcox_Pval, "BH")) %>%
    filter(., !is.nan(Wilcox_Pval))
  
  OY_mean <- {{df}} %>%
    group_by(Age, {{grp}}) %>%
    summarise(., Age_avg = mean({{var}}))
  
  OY_FC <- OY_mean %>%
    group_by({{grp}}) %>%
    group_modify(~ {
      .x %>%
        mutate(., log2FC_OvY = log2(Age_avg/Age_avg[Age == "Young"]))
    }) %>%
    filter(., log2FC_OvY !=0)
  
  stat.all <- bind_cols(OY_stat, OY_FC) %>%
    select(., -matches("Age|Age_Avg")) %>%
    mutate(., Direction = ifelse(log2FC_OvY>0, "Old", "Young"))
  return(stat.all)
}
