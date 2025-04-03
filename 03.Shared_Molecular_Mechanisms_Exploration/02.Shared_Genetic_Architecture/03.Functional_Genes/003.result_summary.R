library(data.table)
library(dplyr)

#trait_anxiety <- c("anxiety2021","anxiety2018")
#trait_allergy <- c("BAD","ADE2015","ADE2021","AR_Diagnose","AR_SelfReport","Asthma2020")
trait_anxiety <- c("anxiety2021")
trait_allergy <- c("BAD","AR_Diagnose","Asthma2020")
#tissue <- c("Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Cortex","WHB","Lung")
tissue <- "Brain_Frontal_Cortex_BA9"

###########procedure: fdr in each trait, intersect#################
result_summary1 <- function(anxiety_dat,allergy_dat,tissue)
{
	anxiety_tmp <- fread(paste0(anxiety_dat,"_",tissue,".smr")) %>% as.data.frame() %>% select(probeID,ProbeChr,Gene,Probe_bp,b_SMR,se_SMR,p_SMR,p_HEIDI,nsnp_HEIDI)
	allergy_tmp <- fread(paste0(allergy_dat,"_",tissue,".smr")) %>% as.data.frame() %>% select(probeID,ProbeChr,Gene,Probe_bp,b_SMR,se_SMR,p_SMR,p_HEIDI,nsnp_HEIDI)
	
	anxiety_tmp$fdr_adj_p <- p.adjust(anxiety_tmp$p_SMR, method = "fdr")
	allergy_tmp$fdr_adj_p <- p.adjust(allergy_tmp$p_SMR, method = "fdr")
	
	anxiety_tmp_sig <- subset(anxiety_tmp,fdr_adj_p<0.05 & p_HEIDI>0.05)
	allergy_tmp_sig <- subset(allergy_tmp,fdr_adj_p<0.05 & p_HEIDI>0.05)
	if(nrow(anxiety_tmp_sig)==0 | nrow(allergy_tmp_sig)==0)
	{
		return("can not be intersected")
	}
	else
	{
		allergy_tmp_sig_gene <- data.frame(allergy_tmp_sig$Gene)
		colnames(allergy_tmp_sig_gene) <- "Gene"
		merge_gene <- merge(anxiety_tmp_sig,allergy_tmp_sig_gene,by="Gene")
	}

	if(nrow(merge_gene) == 0)
	{
		return("can not be intersected")
	}
	else
	{
		merge_gene$tissue <- tissue
		merge_gene$trait1 <- anxiety_dat
		merge_gene$trait2 <- allergy_dat
		merge_gene <- select(merge_gene,trait1,trait2,probeID,ProbeChr,Gene,Probe_bp,b_SMR,se_SMR,p_SMR,fdr_adj_p,p_HEIDI,nsnp_HEIDI,tissue)
		return(merge_gene)
	}
}

########procedure: intersect, fdr########
result_summary2 <- function(anxiety_dat,allergy_dat,tissue)
{
        anxiety_tmp <- fread(paste0(anxiety_dat,"_",tissue,".smr")) %>% as.data.frame() %>% select(probeID,ProbeChr,Gene,Probe_bp,b_SMR,se_SMR,p_SMR,p_HEIDI,nsnp_HEIDI)
        allergy_tmp <- fread(paste0(allergy_dat,"_",tissue,".smr")) %>% as.data.frame() %>% select(probeID,ProbeChr,Gene,Probe_bp,b_SMR,se_SMR,p_SMR,p_HEIDI,nsnp_HEIDI)

        #anxiety_tmp$fdr_adj_p <- p.adjust(anxiety_tmp$p_SMR, method = "fdr")
        #allergy_tmp$fdr_adj_p <- p.adjust(allergy_tmp$p_SMR, method = "fdr")

        #anxiety_tmp_sig <- subset(anxiety_tmp,fdr_adj_p<0.05 & p_HEIDI>0.05)
        #allergy_tmp_sig <- subset(allergy_tmp,fdr_adj_p<0.05 & p_HEIDI>0.05)
        anxiety_tmp_sig <- anxiety_tmp
	allergy_tmp_sig <- allergy_tmp
	if(nrow(anxiety_tmp_sig)==0 | nrow(allergy_tmp_sig)==0)
        {
                return("can not be intersected")
        }
        else
        {
                allergy_tmp_sig_gene <- data.frame(allergy_tmp_sig$Gene)
                colnames(allergy_tmp_sig_gene) <- "Gene"
                merge_gene <- merge(anxiety_tmp_sig,allergy_tmp_sig_gene,by="Gene")
        }

        if(nrow(merge_gene) == 0)
        {
                return("can not be intersected")
        }
        else
        {
                merge_gene$tissue <- tissue
                merge_gene$trait1 <- anxiety_dat
                merge_gene$trait2 <- allergy_dat

		merge_gene$fdr_adj_p <- p.adjust(merge_gene$p_SMR, method = "fdr")
		merge_gene <- subset(merge_gene,fdr_adj_p<0.05 & p_HEIDI>0.05)

		if(nrow(merge_gene) == 0)
		{
			return("can not be intersected")
		}
		else
		{
                merge_gene <- select(merge_gene,trait1,trait2,probeID,ProbeChr,Gene,Probe_bp,b_SMR,se_SMR,p_SMR,fdr_adj_p,p_HEIDI,nsnp_HEIDI,tissue)
                return(merge_gene)
		}
        }
}

########procedure: intersect, BF########
result_summary3 <- function(anxiety_dat,allergy_dat,tissue)
{
        anxiety_tmp <- fread(paste0(anxiety_dat,"_",tissue,".smr")) %>% as.data.frame() %>% select(probeID,ProbeChr,Gene,Probe_bp,b_SMR,se_SMR,p_SMR,p_HEIDI,nsnp_HEIDI)
        allergy_tmp <- fread(paste0(allergy_dat,"_",tissue,".smr")) %>% as.data.frame() %>% select(probeID,ProbeChr,Gene,Probe_bp,b_SMR,se_SMR,p_SMR,p_HEIDI,nsnp_HEIDI)

        #anxiety_tmp$fdr_adj_p <- p.adjust(anxiety_tmp$p_SMR, method = "fdr")
        #allergy_tmp$fdr_adj_p <- p.adjust(allergy_tmp$p_SMR, method = "fdr")

        #anxiety_tmp_sig <- subset(anxiety_tmp,fdr_adj_p<0.05 & p_HEIDI>0.05)
        #allergy_tmp_sig <- subset(allergy_tmp,fdr_adj_p<0.05 & p_HEIDI>0.05)
        anxiety_tmp_sig <- anxiety_tmp
        allergy_tmp_sig <- allergy_tmp
	for(i in 5:9)
	{
		colnames(anxiety_tmp_sig)[i] <- paste0(colnames(anxiety_tmp_sig)[i],"_trait1")
		colnames(allergy_tmp_sig)[i] <- paste0(colnames(allergy_tmp_sig)[i],"_trait2")
	}
        if(nrow(anxiety_tmp_sig)==0 | nrow(allergy_tmp_sig)==0)
        {
                return("can not be intersected")
        }
        else
        {
                merge_gene <- merge(anxiety_tmp_sig,allergy_tmp_sig,by=c("probeID","ProbeChr","Gene","Probe_bp"))
        }

        if(nrow(merge_gene) == 0)
        {
                return("can not be intersected")
        }
        else
        {
                merge_gene$tissue <- tissue
                merge_gene$trait1 <- anxiety_dat
                merge_gene$trait2 <- allergy_dat

                merge_gene$FDR_trait1 <- p.adjust(merge_gene$p_SMR_trait1, method = "fdr")
		merge_gene$FDR_trait2 <- p.adjust(merge_gene$p_SMR_trait2, method = "fdr")

                merge_gene <- subset(merge_gene,p_SMR_trait1<0.05 & p_SMR_trait2<0.05 & p_HEIDI_trait1>0.05 & p_HEIDI_trait2>0.05)

                if(nrow(merge_gene) == 0)
                {
                        return("can not be intersected")
                }
                else
                {
                merge_gene <- select(merge_gene,trait1,trait2,probeID,ProbeChr,Gene,Probe_bp,b_SMR_trait1,se_SMR_trait1,p_SMR_trait1,FDR_trait1,p_HEIDI_trait1,nsnp_HEIDI_trait1,b_SMR_trait2,se_SMR_trait2,p_SMR_trait2,FDR_trait2,p_HEIDI_trait2,nsnp_HEIDI_trait2,tissue)
                return(merge_gene)
                }
        }
}



info <- data.frame()
for(k in 1:length(trait_anxiety))
{
	for(i in 1:length(trait_allergy))
	{
		for(j in 1:length(tissue))
		{
			temp <- result_summary3(trait_anxiety[k],trait_allergy[i],tissue[j])
			if(class(temp) == "character")
			{
				next
			}
			else
			{
				if(nrow(info) == 0)
				{
					info <- temp
				}
				else
				{
					info <- rbind(info,temp)
				}
			}
		}
	}
}
fwrite(info,"info20241211_p005_BFC.txt",sep="\t")
