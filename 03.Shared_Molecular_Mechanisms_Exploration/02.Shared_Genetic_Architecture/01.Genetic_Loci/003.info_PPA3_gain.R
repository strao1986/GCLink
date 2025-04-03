library(data.table)
library(dplyr)

data_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/03.GWAS-PW/"

anxiety_trait <- c("anxiety2018","anxiety2021")
AD_trait <- c("BAD","SelfReportAR","DiagnosedAR","Asthma2020","Asthma2018","ADE2021","ADE2015")

info_PPA3 <- data.frame()
for(i in 1:length(anxiety_trait))
{
  for(j in 1:length(AD_trait))
  {
    filename <- paste0(anxiety_trait[i],"_",AD_trait[j])
    setwd(paste0(data_dir,filename))
    result_SNPs <- fread(paste0(filename,".bfs"))
    result_regions <- fread(paste0(filename,".segbfs"))
    info_temp <- data.frame(phenotype1 = anxiety_trait[i], phenotype2 = AD_trait[j],
                            No.of_SNPs_PPA3_0.8 = length(which(result_SNPs$PPA_3 > 0.8)),
                            No.of_SNPs_PPA3_0.9 = length(which(result_SNPs$PPA_3 > 0.9)),
                            No.of_regions_PPA3_0.8 = length(which(result_regions$PPA_3 > 0.8)),
                            No.of_regions_PPA3_0.9 = length(which(result_regions$PPA_3 > 0.9)))
    if(nrow(info_PPA3) == 0)
    {
      info_PPA3 <- info_temp
    }
    else
    {
      info_PPA3 <- rbind(info_PPA3,info_temp)
    }
  }
}

fwrite(info_PPA3,paste0(data_dir,"info_outline20240909.csv"))


rm(list = ls())
data_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/03.GWAS-PW/"
anxiety_trait <- c("anxiety2021","anxiety2018")
AD_trait <- c("BAD","SelfReportAR","DiagnosedAR","Asthma2020","Asthma2018")
info_PPA3_08_SNP <- data.frame()
info_PPA3_08_loci <- data.frame()
for(i in 1:length(anxiety_trait))
{
  for(j in 1:length(AD_trait))
  {
    filename <- paste0(anxiety_trait[i],"_",AD_trait[j])
    setwd(paste0(data_dir,filename))
    
    result_SNPs <- fread(paste0(filename,".bfs"))
    result_SNPs$Trait1 <- anxiety_trait[i]
    result_SNPs$Trait2 <- AD_trait[j]
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


for(i in 1:length(anxiety_trait))
{
  for(j in 1:length(AD_trait))
  {
    filename <- paste0(anxiety_trait[i],"_",AD_trait[j])
    setwd(paste0(data_dir,filename))
    result_loci <- fread(paste0(filename,".segbfs"))
    result_loci$Trait1 <- anxiety_trait[i]
    result_loci$Trait2 <- AD_trait[j]
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
fwrite(info_PPA3_08_SNP,paste0(data_dir,"info_PPA3_08_SNP20240909.csv"))
fwrite(info_PPA3_08_loci,paste0(data_dir,"info_PPA3_08_loci20240909.csv"))

##################nearest gene annotation#####################
#####generate Annovar input file#####
rm(list = ls())

data_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/03.GWAS-PW/"
setwd(data_dir)
SNP_PPA3_08 <- fread("info_PPA3_08_SNP20240909.csv")

anxiety2021_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/data/clean_data/"
anxiety2021 <- fread(paste0(anxiety2021_dir,"anxiety2021.clean.txt"))

merge_data <- merge(SNP_PPA3_08,select(anxiety2021,pos,chr,rsid,A1,A2),by.x = c("id","pos"),by.y = c("rsid","pos"))
merge_data <- select(merge_data,Trait1,Trait2,chr.y,id,pos,A1,A2,PPA_3)

annovar_input <- select(merge_data,chr.y,pos,A2,A1,id,Trait1,Trait2,PPA_3)
colnames(annovar_input) <- c("chr","pos","A2","A1","rsid","Trait1","Trait2","PPA3")
annovar_input$pos1 <- annovar_input$pos
annovar_input <- select(annovar_input,chr,pos,pos1,A2,A1,rsid,Trait1,Trait2,PPA3)
fwrite(annovar_input,paste0(data_dir,"annovar_input20240909.txt"),sep = "\t",col.names = F)

######gene-level annotation via ANNOVAR in UBUNTU##########
setwd(paste0(data_dir,"Annovar_SNP_annotation/annovar"))
annovar_result <- fread("gene_annotation20240909.variant_function")
annovar_result <- subset(annovar_result,V10 != "Asthma2018")
annovar_result <- select(annovar_result,V8,V3,V4,V7,V6,V2,V9,V10,V11)
colnames(annovar_result) <- c("SNP","Chromosome","BP position","Effect allele","Reference allele","Nearest gene","Trait1","Trait2","PPA3")

fwrite(annovar_result,paste0(data_dir,"gene_annotation20240909.csv"))


data_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/03.GWAS-PW/"
SNP_PPA3_08 <- fread(paste0(data_dir,"info_PPA3_08_SNP20240909.csv"))
#exploration datasets
anxiety2021_DiagnosedAR <- subset(SNP_PPA3_08,Trait1=="anxiety2021" & Trait2=="DiagnosedAR")
anxiety2021_asthma2020 <- subset(SNP_PPA3_08,Trait1=="anxiety2021" & Trait2=="Asthma2020")
anxiety2021_BAD <- subset(SNP_PPA3_08,Trait1=="anxiety2021" & Trait2=="BAD")
#validation datasets
anxiety2021_SelfreportedAR <- subset(SNP_PPA3_08,Trait1=="anxiety2021" & Trait2=="SelfReportAR")
anxiety2018_DiagnosedAR <- subset(SNP_PPA3_08,Trait1=="anxiety2018" & Trait2=="DiagnosedAR")
anxiety2018_asthma2020 <- subset(SNP_PPA3_08,Trait1=="anxiety2018" & Trait2=="Asthma2020")
anxiety2018_BAD <- subset(SNP_PPA3_08,Trait1=="anxiety2018" & Trait2=="BAD")

#AD&AR
intersect(anxiety2021_DiagnosedAR$id,anxiety2021_SelfreportedAR$id) #6
intersect(anxiety2021_DiagnosedAR$id,anxiety2018_DiagnosedAR$id) #9
intersect(intersect(anxiety2021_DiagnosedAR$id,anxiety2021_SelfreportedAR$id),intersect(anxiety2021_DiagnosedAR$id,anxiety2018_DiagnosedAR$id))#2
#AD&asthma
intersect(anxiety2021_asthma2020$id,anxiety2018_asthma2020$id) #15
#AD&BAD
intersect(anxiety2021_BAD$id,anxiety2018_BAD$id) #14


overlap_AR <- merge(subset(SNP_PPA3_08,Trait2=="DiagnosedAR"),subset(SNP_PPA3_08,Trait2=="SelfReportAR"),by = c("id","chr","pos"))
overlap_asthma <- merge(subset(SNP_PPA3_08,Trait2=="Asthma2020"),subset(SNP_PPA3_08,Trait2=="Asthma2018"),by = c("id","chr","pos"))
colnames(overlap_AR) <- c("rsid","chr","pos","Trait1","Trait2","PPA3","Trait1_validation","Trait2_validation","PPA3_validation")
colnames(overlap_asthma) <- c("rsid","chr","pos","Trait1","Trait2","PPA3","Trait1_validation","Trait2_validation","PPA3_validation")
overlap_AR_asthma <- rbind(overlap_AR,overlap_asthma)
overlap_AR_asthma$chr <- as.integer(gsub("chr","",overlap_AR_asthma$chr))

anxiety2021_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/data/clean_data/"
anxiety2021 <- fread(paste0(anxiety2021_dir,"anxiety2021.clean.txt"))
overlap_AR_asthma_full <- merge(overlap_AR_asthma,anxiety2021[,1:5],by=c("chr","pos","rsid"))
overlap_AR_asthma_full <- select(overlap_AR_asthma_full,rsid,chr,pos,A1,A2,everything())
overlap_AR_asthma_full <- arrange(overlap_AR_asthma_full,Trait2)
fwrite(overlap_AR_asthma_full,paste0(data_dir,"info_full20240822.csv"))
