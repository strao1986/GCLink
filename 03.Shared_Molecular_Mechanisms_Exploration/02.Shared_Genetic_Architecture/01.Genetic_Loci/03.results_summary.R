library(data.table)
library(dplyr)

rm(list = ls())

GWASPW_results_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/Github_GCLink/03.Shared_Molecular_Mechanisms_Exploration/02.Shared_Genetic_Architecture/01.Genetic_Loci/" #modify into your path

pheno1 <- c("Anxiety2021")
pheno2 <- c("DiagnosedAR")

#### results summary ####
### SNP ###
info_PPA3_08_SNP <- data.frame()
info_PPA3_08_loci <- data.frame()
for(i in 1:length(pheno1))
{
  for(j in 1:length(pheno2))
  {
    filename <- paste0(pheno1[i],"_",pheno2[j])
    setwd(GWASPW_results_dir)
    
    result_SNPs <- fread(paste0(filename,".bfs"))
    result_SNPs$Trait1 <- pheno1[i]
    result_SNPs$Trait2 <- pheno2[j]
    result_SNPs_PPA3_08 <- subset(result_SNPs,result_SNPs$PPA_3 > 0.8)
    result_SNPs_PPA3_08 <- select(result_SNPs_PPA3_08,Trait1,Trait2,id,chr,pos,PPA_3)
    if(nrow(info_PPA3_08_SNP) == 0)
    {
      info_PPA3_08_SNP <- result_SNPs_PPA3_08
    }
    else
    {
      info_PPA3_08_SNP <- rbind(info_PPA3_08_SNP,result_SNPs_PPA3_08)
    }
  }
}

### loci ###
for(i in 1:length(pheno1))
{
  for(j in 1:length(pheno2))
  {
    filename <- paste0(pheno1[i],"_",pheno2[j])
    setwd(GWASPW_results_dir)
    
    result_loci <- fread(paste0(filename,".segbfs"))
    result_loci$Trait1 <- pheno1[i]
    result_loci$Trait2 <- pheno2[j]
    result_loci_PPA3_08 <- subset(result_loci,result_loci$PPA_3 > 0.8)
    result_loci_PPA3_08 <- select(result_loci_PPA3_08,Trait1,Trait2,chunk,NSNP,chr,st,sp,PPA_3)
    if(nrow(result_loci_PPA3_08) == 0)
    {
      next
    }
    else
    {
      if(nrow(info_PPA3_08_loci) == 0)
      {
        info_PPA3_08_loci <- result_loci_PPA3_08
      }
      else
      {
        info_PPA3_08_loci <- rbind(info_PPA3_08_loci,result_loci_PPA3_08)
      }
    }
  }
}
fwrite(info_PPA3_08_SNP,paste0(GWASPW_results_dir,"info_PPA3_08_SNP.csv"))
fwrite(info_PPA3_08_loci,paste0(GWASPW_results_dir,"info_PPA3_08_loci.csv"))

