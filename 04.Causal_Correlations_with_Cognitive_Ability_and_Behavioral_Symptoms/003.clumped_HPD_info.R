library(data.table)
library(dplyr)
library(TwoSampleMR)

rm(list = ls())
gc()
#################################clumped^exp_SNP information##################################
exp_trait1_path <- 'E:/Projects_documents/data/GWAS_sumstat20241101/EUR/behavior_15/'
exp_trait2_path <- 'E:/Projects_documents/data/GWAS_sumstat20241101/EUR/cognition_8/'

setwd(exp_trait1_path)
exp_trait1 <- str_split_fixed(list.files(pattern = "clean.txt"),"\\.",3)[,1]
exp_trait2 <- "gFactor"
exp_trait <- append(exp_trait1,exp_trait2)
out_trait <- c("Anxiety2021","Anxiety2018","BAD","DiagnosedAR","SelfReportAR","Asthma2020")

info_clumped <- data.frame()
datapath <- 'E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/01.MR/02.5e-6/5.0.Behavior_cognition_to_traits/'
#clump_dir <- "/public/jiangjw/02.allergy_anxiety/data/clump_result/02.p5e6/"

outtran_tmp <- out_trait[1]  #only focus on the number of exposure clumped SNPs, [1] is a casual selection 
for(i in 1:length(exp_trait))
{
  clumped_snp <- fread(paste0(datapath,exp_trait[i],"-",outtran_tmp,"/clump^exp_clean.csv"))
  info_clumped[i,1] <- exp_trait[i]
  info_clumped[i,2] <- nrow(clumped_snp)
}

colnames(info_clumped) <- c("exposure","clumped_snp")
#fwrite(info,paste0(datapath,"info_clumped&MRA_snp.csv"))

##################heterogeneity,pleiotropy,directionality#################
setwd(datapath)
source("MRsummary.func.R")
info_HPD <- data.frame()

exptran <- exp_trait
outtran <- out_trait

for(i in 1:length(exptran))
{
  for(j in 1:length(outtran))
  {
  	if(length(list.files(paste0(datapath,exptran[i],"-",outtran[j]),pattern = "heterogeneity")) == 0)
  	{
  		next
  	}
  	mr_result <- fread(paste0(datapath,exptran[i],"-",outtran[j],"/MRresult.csv"))
  	heterogeneity <- fread(paste0(datapath,exptran[i],"-",outtran[j],"/heterogeneity.csv"))
  	pleiotropy <- fread(paste0(datapath,exptran[i],"-",outtran[j],"/pleiotropy.csv"))
  	directionality <- fread(paste0(datapath,exptran[i],"-",outtran[j],"/directionality_test.csv"))
  	temp <- MR_summary(mr_result,heterogeneity,pleiotropy,directionality)
  	temp1 <- data.frame(temp$summary,Adjusted_pval = temp$summary$`MR pvalue`) ####注意这里的adjusted_pval的p值此时并没有校正，只是为了接下来的校正做准备
  	if(nrow(info_HPD) == 0)
  	{
    		info_HPD <- temp1
  	}
  	else
  	{
    		info_HPD <- rbind(info_HPD,temp1)
  	}


  }
}
#####取iv=1的暴露的结果####
for(i in 1:length(exptran))
{
  for(j in 1:length(outtran))
  {
  	if(length(list.files(paste0(datapath,exptran[i],"-",outtran[j]),pattern = "heterogeneity")) == 0 && length(list.files(paste0(datapath,exptran[i],"-",outtran[j]),pattern = "MRresult.csv")) == 1)
  	{
    		mr_result <- fread(paste0(datapath,exptran[i],"-",outtran[j],"/MRresult.csv"))
    		Wald <- subset(mr_result,method == "Wald ratio")
    		temp2 <- data.frame("exposure"=Wald$exposure,"outcome"=Wald$outcome,"No..of.SNPs.in.MRA"=Wald$nsnp,"Beta"=Wald$b,"Se"=Wald$se,"MR.pvalue"=Wald$pval,"Method"=Wald$method,"heterogeneity"="-","pleiotropy"="-","directionality"="-","OR"=Wald$or,"LCI95"=Wald$or_lci95,"UCI95"=Wald$or_uci95,"Adjusted_pval"=Wald$pval) #这里的adjust_pval也没经过校正
    		info_HPD <- rbind(info_HPD,temp2)
  	}	
  }
}

##############################full information##################################
info_full <- merge(info_clumped,info_HPD,by='exposure')
###############这里校正前一定要留意一下有没有出现weighted median方法与MR egger同时存在的情况###############
#info_full$Adjusted_pval <- p.adjust(info_full$Adjusted_pval,method = "fdr",n = nrow(info_full))   ##到这里统一校正
info_full$beta_se <- paste0(sprintf("%0.3f", info_full$Beta),"(",sprintf("%0.3f", info_full$Se),")")
info_full$or_95CI <- paste0(sprintf("%0.3f", info_full$OR),"(",sprintf("%0.3f", info_full$LCI95),",",sprintf("%0.3f", info_full$UCI95),")")
#beta_se版本
info_full <- select(info_full,exposure,outcome,clumped_snp,No..of.SNPs.in.MRA,MR.pvalue,beta_se,Method,heterogeneity,pleiotropy,directionality,everything())
#OR_95CI版本（二选一）
info_full <- select(info_full,exposure,outcome,clumped_snp,No..of.SNPs.in.MRA,MR.pvalue,or_95CI,Method,heterogeneity,pleiotropy,directionality,everything())

info_full <- arrange(info_full,exposure,outcome)
#info_full$exposure[which(info_full$exposure == "CD4_8")] <- "CD4/CD8"
fwrite(info_full,paste0(datapath,"info_full.csv"))
fwrite(subset(info_full,MR.pvalue<0.05 & (directionality=="TRUE" | directionality=="-")),paste0(datapath,"info_p005_TRUE.csv"))

###############################plot##############################################
for(i in 1:length(exptran))
{
	for(j in 1:length(outtran))
	{
		result_path_plot <- paste0(datapath,as.character(exptran[i]),"-",outtran[j])
  		setwd(result_path_plot)
  		dat_plot <- read.csv("dat_plot.csv")
  		mr_result_plot <- read.csv("MRresult_plot.csv")
  		res_single_plot <- read.csv("res_single_plot.csv")
  		leave_one_out_plot <- read.csv("leave_one_out.csv")

		#scatter plot
		tmp_method_choose <- subset(info_full,exposure==exptran[i] & outcome==outtran[j])
		if(nrow(tmp_method_choose) == 2)
                {
                        pdf(file=paste0("dot_",tmp_method_choose$Method[1],".pdf"))
			print(mr_scatter_plot(subset(mr_result_plot,method == tmp_method_choose$Method[1]),dat_plot))
			dev.off()

			pdf(file=paste0("dot_",tmp_method_choose$Method[2],".pdf"))
                        print(mr_scatter_plot(subset(mr_result_plot,method == tmp_method_choose$Method[2]),dat_plot))
			dev.off()
                }
                else
                {
                        pdf(file='dot.pdf')
			print(mr_scatter_plot(subset(mr_result_plot,method == tmp_method_choose$Method),dat_plot))
			dev.off()
                }
  		pdf(file='forest.pdf', height=15,width=8)#新建一个PDF文件，设置名称、宽高及字体等
  		##森林图
  		print(mr_forest_plot(res_single_plot))
  		dev.off()

  		pdf(file='leave_one_out.pdf', height=15,width=8)#新建一个PDF文件，设置名称、宽高及字体等
  		###逐个剔除检验
  		print(mr_leaveoneout_plot(leave_one_out_plot))
  		dev.off()
  		pdf(file='funnel.pdf', height=6,width=8)#新建一个PDF文件，设置名称、宽高及字体等
  		###漏斗图
  		print(mr_funnel_plot(res_single_plot))
  		dev.off()
	}
}

## supplement ##
info_dir <- "E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/01.MR/02.5e-6/5.0.Behavior_cognition_to_traits/"
info <- fread(paste0(info_dir,"info_full.csv")) %>% as.data.frame()
info_explo <- subset(info,outcome == "Anxiety2021" | outcome == "Asthma2020" | outcome == "BAD" | outcome == "DiagnosedAR") %>% arrange(outcome, exposure, OR)
info_valid <- subset(info,outcome == "Anxiety2018" | outcome == "SelfReportAR") %>% arrange(outcome, exposure, OR)
info_explo_valid <- rbind(info_explo,info_valid)
fwrite(info_explo_valid,paste0(info_dir,"info_explo_valid.csv"))
