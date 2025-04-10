library(data.table)
library(dplyr)
library(tidyverse)

rm(list = ls())
gc()

GWAS_sumstats_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/Github_GCLink/GWAS_Summary_Statistics_Example_Preprocessed/Diseases/" #modify into your path
result_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/Github_GCLink/03.Shared_Molecular_Mechanisms_Exploration/02.Shared_Genetic_Architecture/01.Genetic_Loci/" #modify into your path

pheno1 <- "Anxiety2021"
pheno2 <- "DiagnosedAR"

#### GWAS-pw input making ####
pre_gwas.pw <- function(trait1,trait2,data_dir,result_dir){
  
  for(i in 1:length(trait1)){
    phe1 = fread(paste0(data_dir,trait1[i],".clean.txt")) %>%
      mutate(Z_pheno1 = beta/se) %>%
      mutate(V_pheno1 = se^2) %>%
      rename(SNPID = rsid,
             CHR = chr,
             POS = pos,
      ) %>%
      select(SNPID,CHR,POS,Z_pheno1,V_pheno1)
    for(j in 1:length(trait2)){
      phe2=fread(paste0(data_dir,trait2[j],".clean.txt")) %>% 
        mutate(Z_pheno2 = beta/se) %>%
        mutate(V_pheno2 = se^2) %>%
        rename(SNPID = rsid,    
               CHR = chr,
               POS = pos,
        ) %>%
        select(SNPID,CHR,POS,Z_pheno2,V_pheno2)
      ### 01.extract command SNP(intersection) ###
      gwasp_pw = merge(phe1,phe2,by=c("SNPID","CHR","POS")) #By default the data frames are merged on the columns with names they both have
      ### 02.transform NA into character 'NA' for 4:7 columns of result ###
      ### 03.sort with chr:bpos ###
      gwas_pw_data = gwasp_pw %>%
        group_by(CHR) %>%
        arrange(CHR,POS)
        #mutate_at(vars(contains("V_"),contains("Z_")), as.character) # select the col contains "V_" and"Z_", then transform 
      ### 04.in case se == 0, 'NA' instead of NULL(NA not equal to NUll). ###
      # gwas_pw_data[which(is.na(gwas_pw_data$Z_pheno1)), 4] <- as.character('NA')
      # gwas_pw_data[which(is.na(gwas_pw_data$V_pheno1)), 5] <- as.character('NA')
      # gwas_pw_data[which(is.na(gwas_pw_data$Z_pheno2)), 6] <- as.character('NA')
      # gwas_pw_data[which(is.na(gwas_pw_data$V_pheno2)), 7] <- as.character('NA')
      colnames(gwas_pw_data)[4:7] = c(paste0("Z_",trait1[i]),
                                      paste0("V_",trait1[i]),
                                      paste0("Z_",trait2[j]),
                                      paste0("V_",trait2[j]))
      ### 05. export ###
      #dir.create(result_dir, recursive = T)
      gwas_pw_data$CHR <- paste0("chr",gwas_pw_data$CHR)
      fwrite(gwas_pw_data, file = paste0(result_dir,trait1[i],"_",trait2[j],".pw.txt"), sep = "\t")
    }
    
  }
}

pre_gwas.pw(trait1=pheno1,trait2=pheno2,data_dir=GWAS_sumstats_dir,result_dir=result_dir)
