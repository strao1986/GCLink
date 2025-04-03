library(data.table)
library(dplyr)
library(MAGMA.Celltyping)
library(Seurat)
#library(zellkonverter)

#run based on R-4.3.0

rm(list = ls())
gc()

#if(!require("remotes")) install.packages("remotes")
#remotes::install_github("neurogenomics/MAGMA_Celltyping")


#needed:Brain Frontal Cortex BA9, Brain Hippocampus, Whole Blood, Lung and Brain Cortex
##########################################get available ctd data###################################
Sys.setenv(GITHUB_PAT = "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")
ctd <- MAGMA.Celltyping::get_ctd("ctd_allAIBS",storage_dir=storage_dir) #####cortex(Medial Temporal Gyrus)
ctd <- ewceData::ctd() ####cortex/hippocampus
MAGMA.Celltyping::get_ctd("ctd_allKI",storage_dir=storage_dir) #####cortex, hippocampus, hypothalamus and midbrain
ctd <- MAGMA.Celltyping::get_ctd("ctd_allKI",storage_dir=storage_dir) ###frontal cortex
ctd_allKI <- readRDS("ctd_allKI.rds")
saveRDS(ctd,"ctd_cortex_MTG.rds")
#########################################generate extra celltype dataset#####################
save_CTD_dir <- "E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/04.MAGMA/MAGMA_celltype_enrichment/single_cell_dataset"
#### lung (lung_ctd.rds generate from lung.h5ad file from public server, h5ad data is available from human cell atlas database) ####
#extract counts matrix and meta data from .h5ad file
#lung_data_h5ad <- readH5AD(exp_dir)
#expr_matrix <- assay(lung_data_h5ad, "X")  ##expression matrix
#cell_annotations <- colData(lung_data_h5ad) ##cell annotation infomation
#gene_annotations <- rowData(lung_data_h5ad) ##gene annotation infoamtion
#ctd_file <- "/public/jiangjw/02.allergy_anxiety/10.MAGMA_Celltyping/sc_dataset/lung/lung_ctd.rds"
#ctd <- list(exp = expr_matrix,cell_annotations = as.data.frame(cell_annotations))
#saveRDS(ctd, ctd_file)

#lung_dir <- "/home/rstao/dr.rao/008.Anxiety_allergy/MAGMA_Celltyping_CTD/"
#lung_list <- readRDS(paste0(lung_dir,"lung_ctd.rds"))
#metadata_dir <- "/home/rstao/dr.rao/008.Anxiety_allergy/MAGMA_Celltyping_CTD/"
#lung_metadata <- read.csv(paste0(metadata_dir,"TissueSensitivity-Lung-10x_cell_type_2020-03-12.csv"))
#add 2 levels cell clustering results
#lung_list$annot$level1class <- lung_metadata$annotated_cell_identity.ontology_label
#lung_list$annot$level2class <- lung_metadata$annotated_cell_identity.text
#lung_list$annot$cell_id <- colnames(lung_list$exp)
#generate new list
#lung_list_modify <- list(exp=lung_list$exp,annot=lung_list$annot)
#correct gene symbol
#ref_dir <- "/home/rstao/dr.rao/008.Anxiety_allergy/MAGMA_Celltyping_CTD/"
#if(!file.exists("HGNC_homologene.rpt")){
#  download.file("http://www.informatics.jax.org/downloads/reports/HGNC_AllianceHomology.rpt", destfile=paste0(ref_dir,"HGNC_AllianceHomology.rpt")) #human
#}
#lung_list_modify$exp = EWCE::fix_bad_mgi_symbols(lung_list_modify$exp,mrk_file_path=paste0(ref_dir,"HGNC_AllianceHomology.rpt"))
#saveRDS(lung_list_modify,paste0(lung_dir,"lung_ctd_CorrectGeneSymbol.rds"))

lung_dir <- "/home/rstao/dr.rao/008.Anxiety_allergy/MAGMA_Celltyping_CTD/"
lung_list <- readRDS(paste0(lung_dir,"lung_ctd_CorrectGeneSymbol.rds"))
#filter non-matched genes
exp_lung_CorrectSymbol_dropped = EWCE::drop_uninformative_genes(exp=lung_list$exp,level2annot = lung_list$annot$level2class,input_species = "hsapiens",output_species = "hsapiens",convert_orths=F)
annotLevels = list(level1class=lung_list$annot$level1class,level2class=lung_list$annot$level2class)

ctd_lung <- EWCE::generate_celltype_data(
  exp = exp_lung_CorrectSymbol_dropped,
  annotLevels = annotLevels,
  groupName = "lung20240717", ##generate file named "ctd_lung20240717.rda" if setting so
  no_cores = 1,
  input_species = "hsapiens",
  output_species = "hsapiens",
  savePath=save_CTD_dir)

#### WHB (WHB_seurat.rds file generated from GSE137864 (barcodes,features,matrix files) ####
#establish seurat object based on barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
#...
#extract counts matrix and meta data
WHB_dir <- "/public/jiangjw/02.allergy_anxiety/10.MAGMA_Celltyping/sc_dataset/WHB/"
WHB_seurat <- readRDS(paste0(WHB_dir,"WHB_seurat.rds"))
exp <- as.matrix(GetAssayData(object = WHB_seurat, slot = "counts"))
annot <- WHB_seurat@meta.data
annot <- list(cell_id=colnames(exp),level1class=as.character(annot$orig.ident),level2class=as.character(annot$orig.ident))
WHB_list <- list(exp = exp, annot = annot)
saveRDS(WHB_list,paste0(WHB_dir,"WHB_list20240718.rds"))

WHB_dir <- "E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/04.MAGMA/MAGMA_celltype_enrichment/single_cell_dataset/"
WHB_list <- readRDS(paste0(WHB_dir,"WHB_list20240718.rds"))
ref_dir <- "E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/04.MAGMA/MAGMA_celltype_enrichment/"
#correct gene symbol
WHB_list$exp = EWCE::fix_bad_mgi_symbols(WHB_list$exp,mrk_file_path=paste0(ref_dir,"HGNC_AllianceHomology.rpt"))
#filter non-matched genes
exp_WHB_CorrectSymbol_dropped = EWCE::drop_uninformative_genes(exp=WHB_list$exp,level2annot = WHB_list$annot$level2class,input_species = "hsapiens",output_species = "hsapiens",convert_orths=F)
#establish CTD
annotLevels = list(level1class=WHB_list$annot$level1class)  ##WHB only have level 1 cell-type annotation, so only support 1 level is OK, otherwise error if level 2 is same to level 1
ctd_lung <- EWCE::generate_celltype_data(
  exp = exp_WHB_CorrectSymbol_dropped,
  annotLevels = annotLevels,
  groupName = "WHB20240718", 
  no_cores = 1,
  input_species = "hsapiens",
  output_species = "hsapiens",
  savePath=save_CTD_dir)

#### nasopharynx (human cell atlas, QC in public server) ####
setwd("E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/04.MAGMA/MAGMA_celltype_enrichment/single_cell_dataset/nasopharynx20240722")
exp_nas <- readRDS("exp_qc.rds")
meta_nas <- readRDS("meta_data_qc.rds")
annot <- list(cell_id=meta_nas$NAME,level1class=as.character(meta_nas$Coarse_Cell_Annotations),level2class=as.character(meta_nas$Detailed_Cell_Annotations))
nas_list <- list(exp = exp_nas, annot = annot)

ref_dir <- "E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/04.MAGMA/MAGMA_celltype_enrichment/"
nas_list$exp = EWCE::fix_bad_mgi_symbols(nas_list$exp,mrk_file_path=paste0(ref_dir,"HGNC_AllianceHomology.rpt"))
exp_nas_CorrectSymbol_dropped = EWCE::drop_uninformative_genes(exp=nas_list$exp,level2annot = nas_list$annot$level2class,input_species = "hsapiens",output_species = "hsapiens",convert_orths=F)
annotLevels = list(level1class=nas_list$annot$level1class,level2class=nas_list$annot$level2class)

ctd_nasopharynx <- EWCE::generate_celltype_data(
  exp = exp_nas_CorrectSymbol_dropped,
  annotLevels = annotLevels,
  groupName = "nasopharynx20240723", ##generate file named "ctd_lung20240717.rda" if setting so
  no_cores = 1,
  input_species = "hsapiens",
  output_species = "hsapiens",
  savePath=save_CTD_dir)

#### airway(human cell atlas) ####
#airway
library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)
airway_dir <- "/public/jiangjw/02.allergy_anxiety/10.MAGMA_Celltyping/sc_dataset/airway/"
setwd(airway_dir)
exp_airway <- read.table(paste0(airway_dir,"Raw_exprMatrix_hg19.tsv"),header=T)
rownames(exp_airway) <- exp_airway[,1]
exp_airway <- exp_airway[,-1]
airway_seurat <- CreateSeuratObject(counts = exp_airway)
airway_seurat[["percent.mt"]] <- PercentageFeatureSet(airway_seurat, pattern = "^MT-")
airway_seurat[["percent.rb"]] <- PercentageFeatureSet(airway_seurat, pattern = "^RP[SL]")

#去除数据中的聚类信息保证画出想要的小提琴图
active_te <- as.factor(rep("0",77969)) #77969 is number of cell
names(active_te) <- colnames(airway_seurat)
airway_seurat@active.ident <- active_te

png("qc_violin20240723.png")
VlnPlot(airway_seurat,
        fill.by = "feature", # "feature", "ident"
        features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),
        ncol = 4, pt.size = 0) +
  theme(plot.title = element_text(size=10))
dev.off()
library(DoubletFinder)
library(tidyverse)
library(patchwork)
airway_seurat <- NormalizeData(airway_seurat)
airway_seurat <- FindVariableFeatures(airway_seurat, selection.method = "vst", nfeatures = 2000)
airway_seurat <- ScaleData(airway_seurat)
airway_seurat <- RunPCA(airway_seurat)
airway_seurat <- RunUMAP(airway_seurat, dims = 1:50)
airway_seurat <- FindNeighbors(airway_seurat, dims = 1:50)
airway_seurat <- FindClusters(airway_seurat, resolution = 0.1)
png("umap_PCA50_res01_20240723.png")
DimPlot(airway_seurat, reduction = "umap")
dev.off()
sweep.res.list_pbmc <- paramSweep(airway_seurat, PCs = 1:50, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
png("PK_identification_PCA50_res01_20240723.png")
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
dev.off()
mpK <- as.numeric(as.vector(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
annotations <- airway_seurat$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075 * nrow(airway_seurat@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
airway_seurat <- doubletFinder(airway_seurat, PCs = 1:50,
                               pN = 0.25, pK = mpK, nExp = nExp_poi,
                               reuse.pANN = FALSE, sct = FALSE)
print(head(airway_seurat@meta.data[,8:9]))

names(airway_seurat@meta.data)[8] <- "Double_cell_score"
names(airway_seurat@meta.data)[9] <- "Double_cell_status"
png("umap_DoubleCells_PCA50_res01_20240723.png")
DimPlot(airway_seurat, reduction = "umap",
        group.by = "Double_cell_status")
dev.off()

png("violin_DoubleCells_PCA50_res01_20240723.png")
VlnPlot(airway_seurat, group.by = "Double_cell_status",
        features = c("nCount_RNA", "nFeature_RNA"),
        pt.size = 0, ncol = 2)
dev.off()
saveRDS(airway_seurat,"airway_seurat.rds")
#filter cells
exp <- as.matrix(GetAssayData(object = airway_seurat, slot = "counts"))
annot <- airway_seurat@meta.data
annot_qc <- subset(annot,nFeature_RNA<5000 & nCount_RNA<20000 &
                     percent.mt<30 & percent.rb<15 & Double_cell_status=="Singlet")
exp_qc <- exp[,which(colnames(exp) %in% rownames(annot_qc))]

#add annotation info
meta_data <- fread("meta.tsv")
meta_data$Id <- gsub("-", ".", meta_data$Id)
meta_data <- meta_data[which(meta_data$Id %in% rownames(annot_qc)),]
annot_qc$Id <- rownames(annot_qc)
merge_data <- merge(annot_qc,meta_data,by="Id")

saveRDS(exp_qc,"Exp_qc.rds")
saveRDS(merge_data,"meta_data.rds")

#establish CTD object
annot <- list(cell_id=meta_data$Id,level1class=as.character(meta_data$CellType))
airway_list <- list(exp = exp_qc, annot = annot)

ref_dir <- "E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/04.MAGMA/MAGMA_celltype_enrichment/"
airway_list$exp = EWCE::fix_bad_mgi_symbols(airway_list$exp,mrk_file_path=paste0(ref_dir,"HGNC_AllianceHomology.rpt"))

exp_airway_CorrectSymbol_dropped = EWCE::drop_uninformative_genes(exp=airway_list$exp,level2annot = airway_list$annot$level1class,input_species = "hsapiens",output_species = "hsapiens",convert_orths=F)
annotLevels = list(level1class=airway_list$annot$level1class)

save_CTD_dir <- "/home/rstao/dr.rao/008.Anxiety_allergy/MAGMA_Celltyping_CTD/airway/"
ctd_nasopharynx <- EWCE::generate_celltype_data(
  exp = exp_airway_CorrectSymbol_dropped,
  annotLevels = annotLevels,
  groupName = "airway20240724", ##generate file named "ctd_lung20240717.rda" if setting so
  no_cores = 1,
  input_species = "hsapiens",
  output_species = "hsapiens",
  savePath=save_CTD_dir)

#### spleen(human cell atlas)####
#流程：下载h5ad数据；用readH5AD函数读入；保存成rds文件；构建CTD；生成rda文件

#### cotex(human cell atlas) ####
setwd("E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/04.MAGMA/MAGMA_celltype_enrichment/single_cell_dataset/")
cortex_scrna <- readRDS("pre-mRNA.rds")

#### atopic dermatitis ####
## atopic dermatitis (7 healthy control samples) ##
# rds文件事先通过服务器构建完毕，主要是提取7个健康对照的表达谱和注释信息
rds_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/04.MAGMA/MAGMA_celltype_enrichment/single_cell_dataset/Atopic_dermatitis_healthy/"
ref_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/04.MAGMA/MAGMA_celltype_enrichment/"

AD_healthy_list <- readRDS(paste0(rds_dir,"AD_list20250225.rds"))
#correct gene symbol (alterative)
AD_healthy_list$exp = EWCE::fix_bad_mgi_symbols(AD_healthy_list$exp,mrk_file_path=paste0(ref_dir,"HGNC_AllianceHomology.rpt"))
#filter non-expressed genes (alterative)
exp_AD_Healthy_CorrectSymbol_dropped = EWCE::drop_uninformative_genes(exp=AD_healthy_list$exp,level2annot = AD_healthy_list$annot$level1class,input_species = "hsapiens",output_species = "hsapiens",convert_orths=F)
#establish CTD
annotLevels = list(level1class=AD_healthy_list$annot$level1class)  
ctd_AD_healthy <- EWCE::generate_celltype_data(
  exp = AD_healthy_list$exp,
  annotLevels = annotLevels,
  groupName = "AD_healthy", 
  no_cores = 1,
  input_species = "hsapiens",
  output_species = "hsapiens",
  savePath=rds_dir)

## atopic dermatitis (7 disease case samples) ##
# rds文件事先通过服务器构建完毕，主要是提取7个AD样本的表达谱和注释信息
rds_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/04.MAGMA/MAGMA_celltype_enrichment/single_cell_dataset/Atopic_dermatitis_disease/"
ref_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/04.MAGMA/MAGMA_celltype_enrichment/"

AD_disease_list <- readRDS(paste0(rds_dir,"AD_disease_list20250225.rds"))
#correct gene symbol (alterative)
AD_healthy_list$exp = EWCE::fix_bad_mgi_symbols(AD_healthy_list$exp,mrk_file_path=paste0(ref_dir,"HGNC_AllianceHomology.rpt"))
#filter non-expressed genes (alterative)
exp_AD_Healthy_CorrectSymbol_dropped = EWCE::drop_uninformative_genes(exp=AD_healthy_list$exp,level2annot = AD_healthy_list$annot$level1class,input_species = "hsapiens",output_species = "hsapiens",convert_orths=F)
#establish CTD
annotLevels = list(level1class=AD_disease_list$annot$level1class)  
ctd_AD_healthy <- EWCE::generate_celltype_data(
  exp = AD_disease_list$exp,
  annotLevels = annotLevels,
  groupName = "AD_disease", 
  no_cores = 1,
  input_species = "hsapiens",
  output_species = "hsapiens",
  savePath=rds_dir)

#### brain frontal cortex -20250318 ####
## 流程参照ADE ##
rds_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/02.MDD_ADs/MAGMA_Celltyping/Brain_frontal_cortex_PMID39475571/"

BFC_HC <- readRDS(paste0(rds_dir,"BFC.rds"))
#establish CTD
annotLevels = list(level1class=AD_healthy_list$annot$level1class)  
ctd_AD_healthy <- EWCE::generate_celltype_data(
  exp = AD_healthy_list$exp,
  annotLevels = annotLevels,
  groupName = "AD_healthy", 
  no_cores = 1,
  input_species = "hsapiens",
  output_species = "hsapiens",
  savePath=rds_dir)

########################################celltype enrichment#############################
trait <- c("ADE2021","MDD") #EUR
#celltype_dataset <- c("ctd_BlueLake2018_FrontalCortexOnly","ctd_lung20240718","ctd_WHB20240718","ctd_spleen","ctd_cortex_MTG","ctd_AD_healthy","ctd_AD_disease") #all available now
celltype_dataset <- c("ctd_AD_disease") #一次仅能运行一个组织，运行完后需要将结果文件迁移至组织对应的文件夹内，才能继续运行下一个组织，循环此过程
ctd_path <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/04.MAGMA/MAGMA_celltype_enrichment/single_cell_dataset/"
#create directory
for(i in 1:length(trait))
{
  for(j in 1:length(celltype_dataset))
  {
    dir.create(paste0("E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/04.MAGMA/MAGMA_celltype_enrichment/MAGMA_Files/",trait[i],".35UP.10DOWN/",celltype_dataset[j]))
  }
}
####copy magma result files to created directory and modify filenames, such as "anxiety2021.35UP.10DOWN.genes.out"(manual)#######

#####perform celltype enrichment#########
#MAGMA.Celltyping::celltype_associations_pipeline function
for(i in 1:length(trait))
{
  for(j in 1:length(celltype_dataset))
  {
    trait_path <- paste0("E:\\Projects_documents\\05.allergic_diseases_mental_disorders\\01.anxiety_ADs\\04.MAGMA\\MAGMA_celltype_enrichment\\MAGMA_Files\\",trait[i],".35UP.10DOWN")
    #result_path <- paste0("E:\\Projects_documents\\05.allergic_diseases_mental_disorders\\01.anxiety_ADs\\04.MAGMA\\MAGMA_celltype_enrichment\\MAGMA_Files\\",trait[i],".35UP.10DOWN\\",celltype_dataset[j])
    #ctd <- readRDS(paste0(ctd_path,celltype_dataset[j],".rds"))
    ctd_load <- load(paste0(ctd_path,celltype_dataset[j],".rda"))  ##赋予的是ctd变量
    ctd <- MAGMA.Celltyping::prepare_quantile_groups(ctd = ctd,input_species="hsapiens",
                                                     output_species="hsapiens")
    ################Linear mode celltype enrichment###############
    tryCatch({
      MAGMA.Celltyping::calculate_celltype_associations(
        ctd = ctd,
        #ctd_levels = 1,
        magma_dir = trait_path,
        EnrichmentMode = "Linear",
        ctd_species = "hsapiens",
        force_new = TRUE
      )
    }, error = function(e) {
      message("Error in Linear mode for ", trait[i], ": ", e$message)
    })
    ################top 10% mode####################
    tryCatch({
      MAGMA.Celltyping::calculate_celltype_associations(
        ctd = ctd,
        magma_dir = trait_path,
        EnrichmentMode = "Top 10%",
        ctd_species = "hsapiens",
        force_new = TRUE
      )
    }, error = function(e) {
      message("Error in Top 10% mode for ", trait[i], ": ", e$message)
    })
    #tryCatch({
    #  MAGMA_results <- MAGMA.Celltyping::celltype_associations_pipeline(
    #    magma_dirs = trait_path,
    #    ctd = ctd,
    #    ctd_species = "hsapiens", 
    #    ctd_levels = 1,
    #    ctd_name = "test", ###must be " ", otherwise cannot be recognized 
    #    run_linear = TRUE, 
    #    run_top10 = TRUE,
    #    save_dir = result_path)
    #  MAGMA.Celltyping::merge_results(MAGMA_results = MAGMA_results,species = "human",
    #                                  method = "fdr",save_dir = result_path)
    #}, error = function(e) {
      # 发生错误时，打印错误信息并继续下一次循环
    #  message("Error in processing ", trait[i], " and ", celltype_dataset[j], ": ", e$message)
      #return(NULL)
    #  print(e)
    #})
  }
}

#######################result summary####################
#create directory named relative CTD in E:/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/04.MAGMA/MAGMA_celltype_enrichment/MAGMA_Files/trait.35UP.10DOWN 
#move .out.txt file to created CTD directory(all tissues)

#remove ^# rows ↓ 
#run E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/04.MAGMA/MAGMA_celltype_enrichment/MAGMA_Files/001.result_process.sh

#move .out.clean.txt to E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/04.MAGMA/MAGMA_celltype_enrichment/result/

#summary results
#trait <- c("Anxiety2021","Anxiety2018","BAD","SelfReportAR","DiagnosedAR","Asthma2020","Asthma2018","ADE2021","ADE2015")
trait <- c("MDD","ADE2021")
celltype_dataset <- c("ctd_AD_disease")
method <- c("Top10pct","Linear")

data_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/04.MAGMA/MAGMA_celltype_enrichment/result/"

info <- data.frame()
for(i in celltype_dataset)
{
  for(j in trait)
  {
    for(k in method)
    {
        result_dir <- paste0(data_dir,i,"/")
        setwd(result_dir)
        if(length(grep("level2",list.files())) != 0)
        {
          result_tem <- fread(paste0(result_dir,j,".35UP.10DOWN.level2.MainRun.",k,".gsa.out.clean.txt")) %>% as.data.frame()
          result_tem <- select(result_tem,VARIABLE,TYPE,NGENES,BETA,BETA_STD,SE,P)
          result_tem$FDR <- p.adjust(result_tem$P,method = "fdr")
          result_tem$CTD <- i
          result_tem$trait <- j
          result_tem <- subset(result_tem,P<0.05)
          if(nrow(result_tem) == 0)
          {
            next
          }
          else
          {
            if(nrow(info) == 0)
            {
              info <- result_tem
            }
            else
            {
              info <- rbind(info,result_tem)
            }
          } 
        }
        else
        {
          result_tem <- fread(paste0(result_dir,j,".35UP.10DOWN.level1.MainRun.",k,".gsa.out.clean.txt")) %>% as.data.frame()
          result_tem <- select(result_tem,VARIABLE,TYPE,NGENES,BETA,BETA_STD,SE,P)
          result_tem$FDR <- p.adjust(result_tem$P,method = "fdr")
          result_tem$CTD <- i
          result_tem$trait <- j
          result_tem <- subset(result_tem,P<0.05)
          if(nrow(result_tem) == 0)
          {
            next
          }
          else
          {
            if(nrow(info) == 0)
            {
              info <- result_tem
            }
            else
            {
              info <- rbind(info,result_tem)
            }
          }
        }
    }
  }
}
method_annotation <- data.frame(TYPE=c("SET","COVAR"),method=c("Top10%","Linear"))
tissue_annotation <- data.frame(CTD=c("ctd_BlueLake2018_FrontalCortexOnly","ctd_lung20240718","ctd_WHB20240718","ctd_spleen","ctd_cortex_MTG","ctd_AD_healthy","ctd_AD_disease"),
                                tissue=c("Brain frontal cortex (BA9)","Lung","WHB","Spleen","Medial Temporal Gyrus","Skin_HealthySamples","Skin_ADSamples"))
unique(info$VARIABLE)
info <- subset(info,VARIABLE != "Unknown")
unique(info$VARIABLE)
############generate full cell type name by GPT################
# ## Skin_HealthySample ##
# cell_type_abbreviation_trans <- data.frame(abbr=unique(info$VARIABLE),full_name=c("Conventional Dendritic Cell Type 1","Migratory Memory T Cell Type 3","Effector Regulatory T Cell Type 2",
#                                                                                   "Inflammatory Monocyte","Langerhans Cell Type 1","Langerhans Cell Type 2","Langerhans Cell Type 3",
#                                                                                   "Macrophage Type 3","Monocyte-Derived Dendritic Cell Type 3","Migratory Memory T Cell Type 1",
#                                                                                   "Monocyte-Derived Dendritic Cell Type 1","Effector Cytotoxic T Lymphocyte Memory",
#                                                                                   "Exhausted Cytotoxic T Lymphocyte","Conventional Dendritic Cell Type 2",
#                                                                                   "Conventional Dendritic Cell Type 3","Effector Regulatory T Cell Type 1",
#                                                                                   "Innate Lymphoid Cell/Natural Killer Cell","Cycling Mast Cell","Migratory Dendritic Cell",
#                                                                                   "Natural Killer Cell","Plasma Cell","T Effector Cell","Migratory Memory T Cell Type 2","Naive T Cell",
#                                                                                   "Cycling Regulatory T Cell","Cycling Tissue-Resident Memory T Cell","Tissue-Resident Memory T Cell Type 1",
#                                                                                   "Tissue-Resident Memory T Cell Type 2","Tissue-Resident Memory T Cell Type 3"))

## Skin_ADSample ##
cell_type_abbreviation_trans <- data.frame(abbr=unique(info$VARIABLE),full_name=c("Mast Cell (Cycling)","Conventional Dendritic Cell Type 1","Effector Regulatory T Cell Type 2",
                                                                                  "Monocyte-Derived Dendritic Cell Type 2","Plasma Cell","Migratory Memory T Cell Type 3",
                                                                                  "Conventional Dendritic Cell Type 2","Langerhans Cell Type 3","Monocyte-Derived Dendritic Cell Type 1",
                                                                                  "Migratory Memory T Cell Type 1","Tissue-Resident Memory T Cell Type 3","Central Memory Regulatory T Cell",
                                                                                  "Cytotoxic T Lymphocyte (Cycling)","Effector Memory Cytotoxic T Lymphocyte",
                                                                                  "Exhausted Cytotoxic T Lymphocyte","Conventional Dendritic Cell Type 3","Effector Regulatory T Cell Type 1",
                                                                                  "Innate Lymphoid Cell/Natural Killer Cell","Type 2 Innate Lymphoid Cell","Inflammatory Monocyte",
                                                                                  "Langerhans Cell Type 1","Macrophage Type 3","Macrophage Type 4","Migratory Dendritic Cell",
                                                                                  "Central Memory T Cell","Tet Cell","Migratory Memory T Cell Type 2","Naive T Cell",
                                                                                  "Tissue-Resident Memory T Cell Type 1","Tissue-Resident Memory T Cell Type 2"))

## general ##
# cell_type_abbreviation_trans <- data.frame(abbr=unique(info$VARIABLE),full_name=c("Excitatory neuron subtype 4", "Inhibitory neuron subtype 1b", "Inhibitory neuron subtype 8",
#                                                                                   "Excitatory neuron subtype 1", "Excitatory neuron subtype 2", "Excitatory neuron subtype 5b",
#                                                                                   "Excitatory neuron subtype 8", "Inhibitory neuron subtype 1a", "Inhibitory neuron subtype 1c",
#                                                                                   "Inhibitory neuron subtype 3", "Inhibitory neuron subtype 4a", "Inhibitory neuron subtype 4b",
#                                                                                   "Inhibitory neuron subtype 6b", "Inhibitory neuron subtype 7", "Astrocyte", "Microglia",
#                                                                                   "Activated dendritic cell", "B cell", "Blood vessel cell", "Dendritic cell type 2",
#                                                                                   "Dividing natural killer cell", "Lymph vessel cell", "Mast cell", "Monocyte",
#                                                                                   "Muscle cell", "Plasma cell", "Regulatory T cell", "CD4+ T cell", "CD8+ cytotoxic T cell",
#                                                                                   "Dendritic cell type 1", "Dendritic cell dividing monocyte", "Dividing T cell",
#                                                                                   "Plasmacytoid dendritic cell", "Natural killer cell", "Lymphoid-primed multipotent progenitor",
#                                                                                   "Multipotent lymphoid progenitor", "Mature B cell", "Naive B cell", 
#                                                                                   "Natural killer cell CD160 positive", "IgM producing plasma cell", "Conventional CD4+ T cell",
#                                                                                   "Follicular helper CD4+ T cell", "Regulatory CD4+ T cell", "Cytotoxic CD8+ T cell",
#                                                                                   "Gamma-delta CD8+ T cell", "Dendritic cell type 1", "Activated dendritic cell", 
#                                                                                   "Macrophage", "Mucosal-associated invariant T CD8+ cell", "GABAergic neuron", "Glutamatergic neuron",
#                                                                                   "Microglia"))

merge_data <- merge(info,method_annotation,by="TYPE")
merge_data <- merge(merge_data,tissue_annotation,by="CTD")
merge_data_final <- merge(merge_data,cell_type_abbreviation_trans,by.x="VARIABLE",by.y="abbr")
merge_data_final <- select(merge_data_final,trait,tissue,full_name,VARIABLE,method,NGENES,BETA,SE,P,FDR)
colnames(merge_data_final) <- c("Trait","Cell type dataset (Tissue)","Cell type","Cell type abbreviation","Enrichment method","No. of genes","Beta","SE","P","FDR")
result_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/02.MDD_ADs/MAGMA_Celltyping/Skin_ADSamples/"
fwrite(merge_data_final,paste0(result_dir,"info_p005_20250225.csv"))

info_p005 <- fread(paste0(result_dir,"info_p005_20250225.csv"))
info_FDR005 <- subset(info_p005,FDR<0.05)
fwrite(info_FDR005,paste0(result_dir,"info_FDR005_20250225.csv"))

rm(list = ls())
gc()
################conduct overlap############
# anxiety_trait <- c("Anxiety2021","Anxiety2018")
# allergy_trait <- c("BAD","SelfReportAR","DiagnosedAR","Asthma2020","Asthma2018","ADE2021","ADE2015")
anxiety_trait <- c("MDD")
allergy_trait <- c("ADE2021")
# celltype_dataset <- c("Brain frontal cortex (BA9)","Lung","WHB","Spleen","Medial Temporal Gyrus","Skin_HealthySamples","Skin_ADSamples") #这里要修改成组织
celltype_dataset <- c("Skin_ADSamples")

# #######first:overlap, second:FDR adjust#######
# setwd(result_dir)
# p005 <- fread("info_p005_20241106.csv")
# FDR005 <- p005
# colnames(FDR005)[2] <- "Celltypedataset"
# colnames(FDR005) <- gsub(" ","",colnames(FDR005))
# FDR005 <- subset(FDR005,Enrichmentmethod=="Linear")
# FDR005 <- select(FDR005,Trait,Celltypedataset,Celltype,No.ofgenes,Beta,SE,P)
# 
# result_overlap <- data.frame()
# for(i in celltype_dataset)
# {
#   for(j in anxiety_trait)
#   {
#     for(k in allergy_trait)
#     {
#       temp_anxiety <- subset(FDR005,Trait==j & Celltypedataset==i)
#       colnames(temp_anxiety) <- c("Trait1","Cell_type_dataset","Cell_type","No.of_genes_trait1","Beta_trait1","SE_trait1","P_trait1")
#       temp_allergy <- subset(FDR005,Trait==k & Celltypedataset==i)
#       colnames(temp_allergy) <- c("Trait2","Cell_type_dataset","Cell_type","No.of_genes_trait2","Beta_trait2","SE_trait2","P_trait2")
#       if(nrow(temp_anxiety)==0 | nrow(temp_allergy)==0)
#       {
#         next
#       }
#       else
#       {
#         merge_data <- merge(temp_anxiety,temp_allergy,by=c("Cell_type_dataset","Cell_type"))
#         if(nrow(merge_data)==0)
#         {
#           next
#         }
#         else
#         {
#           merge_data$FDR_trait1 <- p.adjust(merge_data$P_trait1,method = "fdr",n=nrow(merge_data))
#           merge_data$FDR_trait2 <- p.adjust(merge_data$P_trait2,method = "fdr",n=nrow(merge_data))
#           if(nrow(result_overlap)==0)
#           {
#             result_overlap <- merge_data
#           }
#           else
#           {
#             result_overlap <- rbind(result_overlap,merge_data)
#           }
#         }
#       }
#     }
#   }
# }
# result_overlap <- select(result_overlap,Cell_type_dataset,Cell_type,Trait1,Trait2,No.of_genes_trait1,Beta_trait1,SE_trait1,P_trait1,FDR_trait1,No.of_genes_trait2,Beta_trait2,SE_trait2,P_trait2,FDR_trait2)
# result_overlap_FDR005 <- subset(result_overlap,FDR_trait1<0.05 & FDR_trait2<0.05)
# fwrite(result_overlap_FDR005,"result_overlap_FDR005_overlap1_adj2_20240910.csv")

########first:FDR adjust, second:overlap############
result_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/02.MDD_ADs/MAGMA_Celltyping/Skin_ADSamples/"
p005 <- fread(paste0(result_dir,"info_p005_20250225.csv"))
FDR005 <- p005
colnames(FDR005)[2] <- "Celltypedataset"
colnames(FDR005) <- gsub(" ","",colnames(FDR005))
FDR005 <- subset(FDR005,Enrichmentmethod=="Linear")
FDR005 <- select(FDR005,Trait,Celltypedataset,Celltype,No.ofgenes,Beta,SE,P,FDR)

result_overlap <- data.frame()
for(i in celltype_dataset)
{
  for(j in anxiety_trait)
  {
    for(k in allergy_trait)
    {
      temp_anxiety <- subset(FDR005,Trait==j & Celltypedataset==i)
      colnames(temp_anxiety) <- c("Trait1","Cell_type_dataset","Cell_type","No.of_genes_trait1","Beta_trait1","SE_trait1","P_trait1","FDR_trait1")
      temp_allergy <- subset(FDR005,Trait==k & Celltypedataset==i)
      colnames(temp_allergy) <- c("Trait2","Cell_type_dataset","Cell_type","No.of_genes_trait2","Beta_trait2","SE_trait2","P_trait2","FDR_trait2")
      if(nrow(temp_anxiety)==0 | nrow(temp_allergy)==0)
      {
        next
      }
      else
      {
        merge_data <- merge(temp_anxiety,temp_allergy,by=c("Cell_type_dataset","Cell_type"))
        if(nrow(merge_data)==0)
        {
          next
        }
        else
        {
          if(nrow(result_overlap)==0)
          {
            result_overlap <- merge_data
          }
          else
          {
            result_overlap <- rbind(result_overlap,merge_data)
          }
        }
      }
    }
  }
}
result_overlap <- select(result_overlap,Cell_type_dataset,Cell_type,Trait1,Trait2,No.of_genes_trait1,Beta_trait1,SE_trait1,P_trait1,FDR_trait1,No.of_genes_trait2,Beta_trait2,SE_trait2,P_trait2,FDR_trait2)
result_overlap_FDR005 <- subset(result_overlap,FDR_trait1<0.05 & FDR_trait2<0.05)
fwrite(result_overlap,paste0(result_dir,"result_overlap_p005_adj1_overlap2_20250225.csv"))
fwrite(result_overlap_FDR005,paste0(result_dir,"result_overlap_FDR005_adj1_overlap2_20250225.csv"))

# #######first:overlap, second:FDR adjust according to the number of cell type#######
# setwd("E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/04.MAGMA/MAGMA_celltype_enrichment")
# p005 <- fread("info_p005_20240724.csv")
# FDR005 <- p005
# colnames(FDR005)[2] <- "Celltypedataset"
# colnames(FDR005) <- gsub(" ","",colnames(FDR005))
# FDR005 <- subset(FDR005,Enrichmentmethod=="Linear")
# FDR005 <- select(FDR005,Trait,Celltypedataset,Celltype,No.ofgenes,Beta,SE,P)
# 
# result_overlap <- data.frame()
# for(i in celltype_dataset)
# {
#   for(j in anxiety_trait)
#   {
#     for(k in allergy_trait)
#     {
#       temp_anxiety <- subset(FDR005,Trait==j & Celltypedataset==i)
#       colnames(temp_anxiety) <- c("Trait1","Cell_type_dataset","Cell_type","No.of_genes_trait1","Beta_trait1","SE_trait1","P_trait1")
#       temp_allergy <- subset(FDR005,Trait==k & Celltypedataset==i)
#       colnames(temp_allergy) <- c("Trait2","Cell_type_dataset","Cell_type","No.of_genes_trait2","Beta_trait2","SE_trait2","P_trait2")
#       if(nrow(temp_anxiety)==0 | nrow(temp_allergy)==0)
#       {
#         next
#       }
#       else
#       {
#         merge_data <- merge(temp_anxiety,temp_allergy,by=c("Cell_type_dataset","Cell_type"))
#         if(nrow(merge_data)==0)
#         {
#           next
#         }
#         else
#         {
#           merge_data$FDR_trait1 <- p.adjust(merge_data$P_trait1,method = "fdr",n=nrow(merge_data))
#           merge_data$FDR_trait2 <- p.adjust(merge_data$P_trait2,method = "fdr",n=nrow(merge_data))
#           if(nrow(result_overlap)==0)
#           {
#             result_overlap <- merge_data
#           }
#           else
#           {
#             result_overlap <- rbind(result_overlap,merge_data)
#           }
#         }
#       }
#     }
#   }
# }
# result_overlap <- select(result_overlap,Cell_type_dataset,Cell_type,Trait1,Trait2,No.of_genes_trait1,Beta_trait1,SE_trait1,P_trait1,FDR_trait1,No.of_genes_trait2,Beta_trait2,SE_trait2,P_trait2,FDR_trait2)
# result_overlap_p005 <- subset(result_overlap,P_trait1<0.05 & P_trait2<0.05)
# 
# celltype_subdataset <- c("ctd_BlueLake2018_FrontalCortexOnly","ctd_lung20240717","ctd_WHB20240718","ctd_nasopharynx20240723","ctd_airway20240724")
# reusult_FDRmerge <- data.frame()
# for(i in 1:length(celltype_subdataset))
# {
#   result_temp <- subset(result_overlap_p005,Cell_type_dataset == celltype_subdataset[i])
#   result_temp$FDR_trait1 <- p.adjust(result_temp$P_trait1,method = "fdr",n=length(unique(result_temp$Cell_type)))
# }
# result_overlap_p005$FDR_trait1 <- p.adjust(result_overlap_p005$P_trait1,method = "fdr",n=)
# fwrite(result_overlap_FDR005,"result_overlap_FDR005_overlap1_adj2_20240816.csv")

#########result summary############
setwd("E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/04.MAGMA/MAGMA_celltype_enrichment")
result <- fread("result_overlap_FDR005_overlap1_adj2_20240910.csv")

celltype_dataset <- c("ctd_BlueLake2018_FrontalCortexOnly","ctd_lung20240717","ctd_WHB20240718","ctd_cortex_MTG","ctd_nasopharynx20240723","ctd_airway20240724")
anxiety_exploration <- "Anxiety2021"
anxiety_validation <- "Anxiety2018"
allergy_exploration <- c("BAD","DiagnosedAR","Asthma2020")
allergy_validation <- c("SelfReportAR","Asthma2018")

result <- result[-grep("ADE",result$Trait2),]
result_exploration <- subset(result,Trait1%in%anxiety_exploration & Trait2%in%allergy_exploration)
result_validation <- subset(result,Trait1%in%anxiety_validation | Trait2%in%allergy_validation)
result_exploration <- select(result_exploration,Trait1,Trait2,Cell_type_dataset,Cell_type,everything())
result_validation <- select(result_validation,Trait1,Trait2,Cell_type_dataset,Cell_type,everything())
colnames(result_exploration) <- c("Trait1","Trait2","Tissue","Cell type","Number of genes (Trait1)","Beta (Trait1)","SE (Trait1)",
                                  "P (Trait1)","FDR (Trait1)","Number of genes (Trait2)","Beta (Trait2)","SE (Trait2)","P (Trait2)","FDR (Trait2)")
colnames(result_validation) <- colnames(result_exploration)
fwrite(result_exploration,"result_exploration20240910.csv")
fwrite(result_validation,"result_validation20240910.csv")
