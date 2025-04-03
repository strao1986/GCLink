library(data.table)
library(dplyr)

traits <- c("Anxiety2017","Anxiety2018","Anxiety2021","BAD","DiagnosedAR","SelfReportAR","Asthma2020","Asthma2018","ADE2021","ADE2015")

result_summary <- function(traits)
{
	tmp <- fread(paste0(traits,"_multitissue_gene_expr.cell_type_results.txt"),header=T) %>% as.data.frame()
	tmp <- arrange(tmp,Coefficient_P_value)
	tmp$fdr_adj_Pvalue <- p.adjust(tmp$Coefficient_P_value,method="fdr")
	tmp$trait <- traits
	tmp <- select(tmp,trait,Name,Coefficient,Coefficient_std_error,Coefficient_P_value,fdr_adj_Pvalue)
	colnames(tmp)[2] <- "tissue"
	return(tmp)
}

info <- data.frame()
for(i in 1:length(traits))
{
	if(i == 1)
	{
		info <- result_summary(traits[i])	
	}
	else
	{
		info <- rbind(info,result_summary(traits[i]))
	}
}
info_p005 <- subset(info,Coefficient_P_value < 0.05)
#fwrite(info,"info.txt",sep="\t")
fwrite(info_p005,"info_p005.txt",sep="\t")
