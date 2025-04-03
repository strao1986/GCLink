library(data.table)
library(dplyr)

tissue_index <- fread("/public/jiangjw/02.anxiety_ADs/reference/LDSC_SEG_annotation_geneset/tissue_index.csv",header=F) %>% as.data.frame()
colnames(tissue_index) <- c("tissue_name","index")
tissue_index$tissue_name <- gsub(" ","_",tissue_index$tissue_name)
print(head(tissue_index))

tissue_dir <- "/public/jiangjw/02.anxiety_ADs/reference/LDSC_SEG_annotation_geneset/Multi_tissue_gene_expr_1000Gv3_ldscores/"

output_dir <- "/public/jiangjw/02.anxiety_ADs/reference/LDSC_SEG_annotation_geneset/"

ldcts_file <- data.frame()
for(i in 1:nrow(tissue_index))
{
	ldcts_file[i,1] <- tissue_index$tissue_name[i]
	ldcts_file[i,2] <- paste0(tissue_dir,"GTEx.",i,".",",",tissue_dir,"GTEx.control.")
}

fwrite(ldcts_file,paste0(output_dir,"Multi_tissue_gene_expr.ldcts"),sep="\t",col.names=F)
