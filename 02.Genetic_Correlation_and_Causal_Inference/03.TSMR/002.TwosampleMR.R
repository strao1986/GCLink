##### 240508#####

rm(list = ls())
gc()
####TowSampleMR包的安装
#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")

library(data.table)
### include packages: data.table, tsmr, raps 
####source('tsmr.func.ad.R')
library(TwoSampleMR)
library(mr.raps)

############################tsmr.func#####################################
mr_pip <- function(exposure, outcome, iv_plink, exp_trait, out_trait, pvalue, respath,samplesize_exp,samplesize_out){
  
  tmp_exp <- fread(exposure) %>% as.data.frame()  ###################读入暴露的clean.txt
  #tmp_exp$F <- (tmp_exp$beta^2)/(tmp_exp$se^2)
  #tmp_exp <- subset(tmp_exp,F>10)
#  tmp_exp$samplesize <- samplesize_exp ###########加入样本量，每行都统一
  tmp_out <- fread(outcome) %>% as.data.frame()  ###################读入结局的clean.txt
#  tmp_out$samplesize <- samplesize_out ###########加入样本量，每行都统一
  tmp_iv <- fread(iv_plink) %>% as.data.frame() ###clump后的结果（暴露），结局不需要做clump
  
  tmp_exp1 <- tmp_exp[which(tmp_exp$p < pvalue), ]    ###########取p-value小于5e-8的SNP
  ### set result path	
  tmp_respath <- paste0(respath, '/', exp_trait, '-', out_trait)   #######创建结果文件夹
  dir.create(tmp_respath, recursive = T)
  ### input exposure data
  exp_dat <- format_data(    ####该函数用于创建单变量MR所需要的暴露的输入文件，注意这里的暴露并未与clump过后的SNP取交集
    type = 'exposure',       ##########这里的所有character都需要根据数据的实际名称进行修改
    tmp_exp1,
    snp_col = "rsid",
    chr_col = "chr",
    pos_col = "pos",
    effect_allele_col = "A1",   #a1对应的是突变后的等位基因
    other_allele_col = "A2",    #a2对应的是参考等位基因
    #eaf_col = "eaf",
    beta_col = "beta",
    se_col = "se",
    pval_col = "p",
    samplesize_col = "N")	  ###经过format_data整理过后的列名会变为"SNP"
  
  ### IV threshold: p-value  ; r2 = 0.1, kb = 1000 
  #	exp_iv <- clump_data(exp_dat, clump_r2 = 0.1, clump_kb = 1000, pop = 'EAS')
  iv <- as.data.frame(tmp_iv$SNP)    ####提取clump后的SNP
  colnames(iv) <- 'SNP'   
  
  exp_iv <- merge(iv, exp_dat, by = "SNP")   ###将暴露组的输入文件的SNP与clump后的SNP取交集，保留交集的行
  exp_iv$exposure <- exp_trait   ######exp_dat中的exposure皆为'exposure'，故在此做替换
  
  ### write down IV information ###
  write.csv(exp_iv, paste0(tmp_respath, '/', "clump^exp_clean.csv"), row.names = F)  ####用于单变量MR分析的暴露的IVs
  
  ### get IVs in  outcome data
  merge <- merge(exp_iv, tmp_out, by.x = "SNP", by.y = "rsid")#取clump后的暴露iv与clean的结局iv的交集
  
  if (nrow(merge) != 0) {
    out_dat <- format_data(   #用于创建单变量MR所需要的结局的输入文件
      type = "outcome",        
      merge,
      snp_col = "SNP",       #根据的是merge的SNP列名
      chr_col = "chr",
      pos_col = "pos",
      effect_allele_col = "A1",  #突变后的的等位
      other_allele_col = "A2",   #突变前的等位
      eaf_col = "eaf",
      beta_col = "beta",
      se_col = "se",
      pval_col = "p",
      samplesize_col = "N")	##该步骤所得到的out_dat即为与暴露的IVs取交集过后的结局的IVs,因为已经merge过了
    
    out_dat$outcome <- out_trait  #修改'outcome'
    dat_raw <- harmonise_data(exp_iv, out_dat, action = 2)   ###注意该函数用于归一化暴露与结局的IV，归一化的对象是暴露与结局的IVs中的a1与a2，因此若数据中不包含a1与a2，则不需要进行归一化，归一化后会去除部分SNP
    #一般我们推荐使用默认值action=2即可，当然也可以使用action=3，这时候就表示去除所有存在回文结构的SNP
    dat <- as.data.frame(dat_raw[which(dat_raw$mr_keep == "TRUE"),])   ##########dat文件真正用于后续的MR分析
    
    directionality_test <- directionality_test(dat)  ##########方向性检验，输入归一化之后的数据，判断因果方向的准确性#######################
    
    #只有mr_keep是TRUE的SNP才是真正用于计算MR结果的IV。如果mr_keep是FALSE的话，那就说明这个SNP在计算MR结果时会被剔除
    if (length(rownames(dat)) != 0){
      ### perform MR analysis
      mr_results <- mr(dat, method_list = c('mr_wald_ratio', 'mr_egger_regression',
                                            'mr_weighted_median', 'mr_ivw', 'mr_simple_mode',
                                            'mr_weighted_mode', 'mr_ivw_fe'))    #选定方法进行单变量MR分析
      write.csv(mr_results, paste0(tmp_respath, "/", "MRresult_plot.csv"), row.names = F)
      write.csv(dat, paste0(tmp_respath, "/", "dat_plot.csv"), row.names = F)
      
      mr_raps <- as.data.frame(mr.raps( dat$beta.exposure, dat$beta.outcome,  ####mr.raps函数的作用就在于使用"Robust Adjusted Profile Score"方法进行MR分析，其是一种当存在较多弱工具变量时，使MR分析结果无偏差的方法，不过raps结果一般都太过显著了，反而不好解释结果
                                        dat$se.exposure, dat$se.outcome))		
      
      ### add cols so that we can conbine mr & raps result
      mr_results$naive.se <- "NA"
      mr_results$chi.sq.test <- "NA"
      
      mr_raps$outcome <- out_trait
      mr_raps$exposure <- exp_trait
      mr_raps$id.exposure <- dat[1,which(colnames(dat) == "id.exposure")]
      mr_raps$id.outcome <- dat[1,which(colnames(dat) == "id.outcome")]
      mr_raps$method <- "Robust Adjusted Profile Score"
      mr_raps$nsnp <- length(which(dat$mr_keep) == "TRUE")
      
      ### correct col locations by mr_results
      mr_raps <- mr_raps[,c(8, 9, 6, 7, 10, 11, 1, 2, 3, 4, 5)]  #####列排序
      colnames(mr_raps) <- c("id.exposure", "id.outcome",
                             "outcome", "exposure", "method",
                             "nsnp", "b", "se", "pval",
                             "naive.se", "chi.sq.test")
      mr_results <- rbind(mr_results, mr_raps)   ######将"Robust Adjusted Profile Score"方法的结果与先前其他方法进行分析的结果合并
      ### 添加OR/CI95%	，如果是二分类性状的话，需要计算OR，如果是连续型的话就用beta即可
      mr_results <- generate_odds_ratios(mr_results)
      #########"naive.se"与"chi.sq.test"表示raps的结果，"lo_ci"与"up_ci"表示beta值的置信区间,后面就是or值与其置信区间
      
      ### single snp analysis
      res_single <- mr_singlesnp(dat)     ####森林图
      write.csv(res_single, paste0(tmp_respath, "/", "res_single_plot.csv"), row.names = F)
      ### output result
      write.csv(out_dat, paste0(tmp_respath, "/", "outcome.csv"), row.names = F)
      write.csv(dat_raw, paste0(tmp_respath, "/", "harmonize.csv"), row.names = F)
      write.csv(mr_results, paste0(tmp_respath, "/", "MRresult.csv"), row.names = F)
      ##########################
      ### sensitive analysis
      if (length(which(dat$mr_keep) == "TRUE") > 1){
        heter_dat <- mr_heterogeneity(dat)   ###异质性分析
        pleio_dat <- mr_pleiotropy_test(dat) ###水平多效性分析
        leave_one_dat <- mr_leaveoneout(dat) ###leave one out test，目的是看有没有对MR分析结果影响特别大的工具变量
        
        write.csv(heter_dat, paste0(tmp_respath, "/", "heterogeneity.csv"), row.names = F)
        write.csv(pleio_dat, paste0(tmp_respath, "/", "pleiotropy.csv"), row.names = F)
        write.csv(leave_one_dat, paste0(tmp_respath, "/", "leave_one_out.csv"), row.names = F)
        write.csv(res_single, paste0(tmp_respath, "/", "single_snp.result.csv"), row.names = F)
        write.csv(directionality_test, paste0(tmp_respath, "/", "directionality_test.csv"), row.names = F)
      }
    }
  }						
}
#########################################################################
exp_traits_anxiety <- c("BAD")
out_traits_IL <- c("anxiety2018","anxiety2021")  #结局变量
#outtran <- 'D:/Projects_documents/05.allergic_diseases&mental_disorders/01.allergic_diseases&anxiety/data/clean_data/anxiety2021.clean.txt'  


### iv_path
clump_F <- 'E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/data/clean_data/clump_result/02.p5e6/' #clump结果目录

exp_anxiety_path <- 'E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/data/clean_data/'
out_data_path <- 'E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/data/clean_data//'

result_path <- 'E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/01.MR/02.5e-6/4.0.BAD/BAD_anxiety/' #结果存放目录

### TSMR

for(i in 1:length(exp_traits_anxiety))
{
	for(j in 1:length(out_traits_IL))
	{
		        tmp_exp <- paste0(exp_anxiety_path,exp_traits_anxiety[i],".clean.txt")
        		tmp_out <- paste0(out_data_path,out_traits_IL[j],".clean.txt")
			tmp_iv <- paste0(clump_F,exp_traits_anxiety[i],".LD_result.clumped")
			mr_pip(exposure = tmp_exp,    
         			outcome = tmp_out,  
         			iv_plink = tmp_iv,     ###clump结果
         			exp_trait = exp_traits_anxiety[i],
         			out_trait = out_traits_IL[j],
         			pvalue = 1,
         			respath = result_path,
         			#samplesize_exp = tmp_samplesize_exp,
         			#samplesize_out = tmp_samplesize_out
        			 )
	
	}
}



#for(k in 1:length(exp_traits_allergy))
#{
#        for(l in 1:length(out_traits_IL))
#        {
#                        tmp_exp <- paste0(exp_allergy_path,exp_traits_allergy[k],".clean.txt")
#                        tmp_out <- paste0(out_data_path,out_traits_IL[l],".clean.txt")
#                        tmp_iv <- paste0(clump_F,exp_traits_allergy[k],".LD_result.clumped")
#                        mr_pip(exposure = tmp_exp,    
#                                outcome = tmp_out,    
#                                iv_plink = tmp_iv,     ###clump结果
#                                exp_trait = exp_traits_allergy[k],
#                                out_trait = out_traits_IL[l],
#                                pvalue = 1,
#                                respath = result_path,
#                                #samplesize_exp = tmp_samplesize_exp,
#                                #samplesize_out = tmp_samplesize_out
#         			)
#        
#        }
#}

