library(data.table)
library(dplyr)

rm(list = ls())

setwd("E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/00.UKB_allergy_anxiety_epidemiology/001.anxiety_allergy")

#######################################clinical information summary##################################
ukb_data <- fread("UKB_data_20240417.clean.txt") %>% as.data.frame
tem <- vector("numeric",0)
for(i in 2:ncol(ukb_data))
{
  if(strsplit(strsplit(colnames(ukb_data)[i],"-")[[1]][2],"\\.")[[1]][1] == "0")
  {
    tem <- append(tem,i)
  }
}
ukb_data_clean <- ukb_data[,c(1,tem)]
colnames(ukb_data_clean)[2:3] <- c("sex","birth_year")
ukb_data_clean$age <- 2017-ukb_data_clean$birth_year
ukb_data_clean <- select(ukb_data_clean,eid,sex,age,everything())

setwd("E:/move/disk_document/Projects_documents/02.UKB_GWAS/traits_data")
anxiety_GP <- fread("nerves_anxiety_tension_depression_GP.txt")
anxiety_psychiatry <- fread("nerves_anxiety_tension_depression_psychiatrist.txt")
anxiety_professional <- fread("Anxiety_nerves_generalizedanxietydisorder_professional.txt")

######################################extracting atopic dermatitis phenotype#############################
#case-L20 control-remain samples excluding L
ADE <- ukb_data_clean[,c(1:3,grep("41270",colnames(ukb_data_clean)))]
ADE$allergy_status <- -9
ADE$skin_status <- 0
for(i in 4:216)
{
  ADE$allergy_status[grep("L20",ADE[,i])] <- 1
  ADE$skin_status[grep("L",ADE[,i])] <- 1
}
ADE$allergy_status[which(ADE$skin_status == 0)] <- 0
ADE <- ADE[-which(ADE$allergy_status == -9),]
ADE <- ADE[,c(1:3,217)]
fwrite(ADE,"ADE.txt",sep = "\t")

#case-1452 control-remain samples excluding categery "dermatology"
dermatitis <- ukb_data_clean[,c(1:3,grep("20002",colnames(ukb_data_clean)))]
dermatitis$allergy_status <- -9
dermatitis$skin_status <- 0
for(i in 4:37)
{
  dermatitis$allergy_status[which(dermatitis[,i] == 1452)] <- 1
  dermatitis$skin_status[which(dermatitis[,i] == 1452 | dermatitis[,i] == 1453 | dermatitis[,i] == 1454 | 
                                dermatitis[,i] == 1455 | dermatitis[,i] == 1548 | dermatitis[,i] == 1549 | 
                                dermatitis[,i] == 1550 | dermatitis[,i] == 1625 | dermatitis[,i] == 1660 | 
                                dermatitis[,i] == 1661 | dermatitis[,i] == 1667 | dermatitis[,i] == 1680)] <- 1
}
dermatitis$allergy_status[which(dermatitis$skin_status == 0)] <- 0
dermatitis <- dermatitis[-which(dermatitis$allergy_status == -9),]
dermatitis <- dermatitis[,c(1:3,38)]
fwrite(dermatitis,"dermatitis.txt",sep = "\t")
#############################################chi-square test#######################################
#three anxiety phenotypes: nerves_anxiety_tension_depression_GP, nerves_anxiety_tension_depression_psychiatrist, Anxiety_nerves_generalizedanxietydisorder_professional
#three allergic diseases phenotypes: dermatitis

anxiety_pheno <- c("Anxiety_nerves_generalizedanxietydisorder_professional_EUR")
allergic_pheno <- c("AR_doctor_EUR","asthma_doctor_EUR","dermatitis_EUR","BAD_EUR")

###############attention!!!, compare_status=1 = compare_status=2###########
chi_square_test <- function(anxiety_phenotype,allergy_phenotype,compare_status)
{
  chi_square_resultsummary <- data.frame()
  temp_anxiety <- fread(paste0("E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/00.UKB_allergy_anxiety_epidemiology/traits_epidemiology_UKB/",anxiety_phenotype,".txt"))
  temp_allergy <- fread(paste0("E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/00.UKB_allergy_anxiety_epidemiology/traits_epidemiology_UKB/",allergy_phenotype,".txt"))
  if(compare_status == 1)  #compare_status == 1 means compare the prevalence of allergic diseases between anxiety group and control group
  {
    chi_square_table <- matrix(c(length(intersect(temp_anxiety$eid[which(temp_anxiety$anxiety_status == 1)],temp_allergy$eid[which(temp_allergy$allergy_status == 1)])),
                                 length(intersect(temp_anxiety$eid[which(temp_anxiety$anxiety_status == 0)],temp_allergy$eid[which(temp_allergy$allergy_status == 1)])),
                                 length(intersect(temp_anxiety$eid[which(temp_anxiety$anxiety_status == 1)],temp_allergy$eid[which(temp_allergy$allergy_status == 0)])),
                                 length(intersect(temp_anxiety$eid[which(temp_anxiety$anxiety_status == 0)],temp_allergy$eid[which(temp_allergy$allergy_status == 0)]))),nrow = 2,ncol = 2)
  }
  
  if(compare_status == 2)  #compare_status == 2 means compare the prevalence of anxiety between allergic diseases group and control group
  {
    chi_square_table <- matrix(c(length(intersect(temp_anxiety$eid[which(temp_anxiety$anxiety_status == 1)],temp_allergy$eid[which(temp_allergy$allergy_status == 1)])),
                                 length(intersect(temp_anxiety$eid[which(temp_anxiety$anxiety_status == 1)],temp_allergy$eid[which(temp_allergy$allergy_status == 0)])),
                                 length(intersect(temp_anxiety$eid[which(temp_anxiety$anxiety_status == 0)],temp_allergy$eid[which(temp_allergy$allergy_status == 1)])),
                                 length(intersect(temp_anxiety$eid[which(temp_anxiety$anxiety_status == 0)],temp_allergy$eid[which(temp_allergy$allergy_status == 0)]))),nrow = 2,ncol = 2)
  }
  
  chi_square_results <- chisq.test(chi_square_table,correct = F)
  chi_square_resultsummary[1,1] <- anxiety_phenotype
  chi_square_resultsummary[1,2] <- unname(table(temp_anxiety$anxiety_status)[2])
  chi_square_resultsummary[1,3] <- unname(table(temp_anxiety$anxiety_status)[1])
  chi_square_resultsummary[1,4] <- allergy_phenotype
  chi_square_resultsummary[1,5] <- unname(table(temp_allergy$allergy_status)[2])
  chi_square_resultsummary[1,6] <- unname(table(temp_allergy$allergy_status)[1])
  chi_square_resultsummary[1,7] <- unname(chi_square_results$statistic)
  chi_square_resultsummary[1,8] <- unname(chi_square_results$parameter)
  chi_square_resultsummary[1,9] <- unname(chi_square_results$p.value)
  return(chi_square_resultsummary)
}

####compare status = 1 
chi_square_resultsummary_sum <- data.frame()
for(i in anxiety_pheno)
{
  for(j in allergic_pheno)
  {
    if(nrow(chi_square_resultsummary_sum) == 0)
    {
      chi_square_resultsummary_sum <- chi_square_test(anxiety_phenotype = i,allergy_phenotype = j,compare_status = 1)
    }
    else
    {
      chi_square_resultsummary_sum <- rbind(chi_square_resultsummary_sum,chi_square_test(anxiety_phenotype = i,allergy_phenotype = j,compare_status = 1))
    }
  }
}
colnames(chi_square_resultsummary_sum) <- c("anxiety_phenotype","anxiety_ncase","anxiety_ncontrol","allergy_phenotype","allergy_ncase","allergy_ncontrol","X-squared","df","p_value")
fwrite(chi_square_resultsummary_sum,"chisquare_result_summary_EUR20240819.txt",sep = "\t")

####compare status = 2 
chi_square_resultsummary_sum <- data.frame()
for(i in anxiety_pheno)
{
  for(j in allergic_pheno)
  {
    if(nrow(chi_square_resultsummary_sum) == 0)
    {
      chi_square_resultsummary_sum <- chi_square_test(anxiety_phenotype = i,allergy_phenotype = j,compare_status = 2)
    }
    else
    {
      chi_square_resultsummary_sum <- rbind(chi_square_resultsummary_sum,chi_square_test(anxiety_phenotype = i,allergy_phenotype = j,compare_status = 1))
    }
  }
}
colnames(chi_square_resultsummary_sum) <- c("anxiety_phenotype","anxiety_ncase","anxiety_ncontrol","allergy_phenotype","allergy_ncase","allergy_ncontrol","X-squared","df","p_value")
fwrite(chi_square_resultsummary_sum,"chisquare_result_summary.txt",sep = "\t")

##################################################logistic regression###############################################
#three anxiety phenotypes: anxiety_GP, anxiety_professional, anxiety_psychiatry
#three allergic diseases phenotypes: dermatitis

anxiety_pheno <- c("Anxiety_nerves_generalizedanxietydisorder_professional_EUR")
allergic_pheno <- c("AR_doctor_EUR","asthma_doctor_EUR","dermatitis_EUR","BAD_EUR")

logistic_regression <- function(anxiety_phenotype,allergy_phenotype,response_status)
{
  logistic_regression_resultsummary <- data.frame()
  temp_anxiety <- fread(paste0("E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/00.UKB_allergy_anxiety_epidemiology/traits_epidemiology_UKB/",anxiety_phenotype,".txt"))
  temp_allergy <- fread(paste0("E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/00.UKB_allergy_anxiety_epidemiology/traits_epidemiology_UKB/",allergy_phenotype,".txt"))
  logistic_data <- merge(x=temp_anxiety,y=temp_allergy,by=c("eid","sex","age"))
  if(response_status == 1)  #response_status == 1 means regarding anxiety as response variable
  {
    logistic_model <- glm(anxiety_status ~ sex + age + allergy_status, data = logistic_data, family = binomial)
  }
  
  if(response_status == 2)  #response_status == 2 means regarding allergic diseases as response variable
  {
    logistic_model <- glm(allergy_status ~ sex + age + anxiety_status, data = logistic_data, family = binomial)
  }
  OR_CI <- data.frame(exp(cbind(OR=coef(logistic_model),confint(logistic_model,level = 0.95))))
  info <- data.frame(coef(summary(logistic_model)))
  colnames(OR_CI) <- c("OR","95%_Lower","95%_Upper")
  logistic_regression_resultsummary[1,1] <- anxiety_phenotype
  logistic_regression_resultsummary[1,2] <- allergy_phenotype
  logistic_regression_resultsummary[1,3] <- OR_CI$OR[4]
  logistic_regression_resultsummary[1,4] <- OR_CI$`95%_Lower`[4]
  logistic_regression_resultsummary[1,5] <- OR_CI$`95%_Upper`[4]
  logistic_regression_resultsummary[1,6] <- info$z.value[4]
  logistic_regression_resultsummary[1,7] <- info$Pr...z..[4]
  return(logistic_regression_resultsummary)
}

####response status = 1 
logistic_responAnxiety_resultsummary_sum <- data.frame()
for(i in anxiety_pheno)
{
  for(j in allergic_pheno)
  {
    if(nrow(logistic_responAnxiety_resultsummary_sum) == 0)
    {
      logistic_responAnxiety_resultsummary_sum <- logistic_regression(anxiety_phenotype = i,allergy_phenotype = j,response_status = 1)
    }
    else
    {
      logistic_responAnxiety_resultsummary_sum <- rbind(logistic_responAnxiety_resultsummary_sum,logistic_regression(anxiety_phenotype = i,allergy_phenotype = j,response_status = 1))
    }
  }
}
colnames(logistic_responAnxiety_resultsummary_sum) <- c("anxiety_phenotype","allergy_phenotype","OR(allergy)","95% Lower","95% Upper","z value","Pr(>|z|)")
fwrite(logistic_responAnxiety_resultsummary_sum,"logistic_responAnxiety_summary_EUR20240820.txt",sep = "\t")

####response status = 2 
logistic_responAllergy_resultsummary_sum <- data.frame()
for(i in anxiety_pheno)
{
  for(j in allergic_pheno)
  {
    if(nrow(logistic_responAllergy_resultsummary_sum) == 0)
    {
      logistic_responAllergy_resultsummary_sum <- logistic_regression(anxiety_phenotype = i,allergy_phenotype = j,response_status = 2)
    }
    else
    {
      logistic_responAllergy_resultsummary_sum <- rbind(logistic_responAllergy_resultsummary_sum,logistic_regression(anxiety_phenotype = i,allergy_phenotype = j,response_status = 2))
    }
  }
}
colnames(logistic_responAllergy_resultsummary_sum) <- c("anxiety_phenotype","allergy_phenotype","OR(anxiety)","95% Lower","95% Upper","z value","Pr(>|z|)")
fwrite(logistic_responAllergy_resultsummary_sum,"logistic_responAllergy_summary.txt",sep = "\t")

##########################################anxiety prevalence chi-square visualization: column plot################################
library(ggplot2)
library(grid)
library(ggsignif)
setwd("E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/00.UKB_allergy_anxiety_epidemiology/traits_epidemiology_UKB")
anxiety_professional <- fread("Anxiety_nerves_generalizedanxietydisorder_professional.txt")
AR <- fread("AR_doctor_EUR.txt")
asthma <- fread("asthma_doctor_EUR.txt")
dermatitis <- fread("dermatitis_EUR.txt")
BAD <- fread("BAD_EUR.txt")

anxiety_prevalence_calculation <- function(anxiety_pheno,allergy_pheno)
{
  anxiety_prevalance <- data.frame()
  anxiety_temp <- fread(paste0("E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/00.UKB_allergy_anxiety_epidemiology/traits_epidemiology_UKB/",anxiety_pheno,".txt"))
  allergy_temp <- fread(paste0("E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/00.UKB_allergy_anxiety_epidemiology/traits_epidemiology_UKB/",allergy_pheno,".txt"))
  anxiety_prevalance[1,1] <- length(intersect(anxiety_temp$eid[which(anxiety_temp$anxiety_status == 1)],allergy_temp$eid[which(allergy_temp$allergy_status == 1)]))/length(which(allergy_temp$allergy_status == 1))
  anxiety_prevalance[1,2] <- length(intersect(anxiety_temp$eid[which(anxiety_temp$anxiety_status == 1)],allergy_temp$eid[which(allergy_temp$allergy_status == 0)]))/length(which(allergy_temp$allergy_status == 0))
  colnames(anxiety_prevalance) <- c("case_prevalence","control_prevalence")
  return(anxiety_prevalance)
}

traits <- c("AR_doctor_EUR","asthma_doctor_EUR","dermatitis_EUR","BAD_EUR")
anxiety_prevalance_summary <- data.frame()
for(i in 1:length(traits))
{
  if(nrow(anxiety_prevalance_summary) == 0)
  {
    anxiety_prevalance_summary <- anxiety_prevalence_calculation("Anxiety_nerves_generalizedanxietydisorder_professional_EUR",traits[i])
  }
  else
  {
    anxiety_prevalance_summary <- rbind(anxiety_prevalance_summary,anxiety_prevalence_calculation("Anxiety_nerves_generalizedanxietydisorder_professional_EUR",traits[i]))
  }
}

plot_data <- cbind(data.frame(traits),anxiety_prevalance_summary)
plot_data <- data.frame(traits=c(rep("AR",2),rep("Asthma",2),rep("ADE",2),rep("BAD",2)),
                        status=c("Allergic disease case","Non-allergy control","Allergic disease case","Non-allergy control","Allergic disease case","Non-allergy control","Allergic disease case","Non-allergy control"),
                        anxiety_prevalence=c(sprintf("%0.3f",anxiety_prevalance_summary[1,1]),sprintf("%0.3f",anxiety_prevalance_summary[1,2]),
                                             sprintf("%0.3f",anxiety_prevalance_summary[2,1]),sprintf("%0.3f",anxiety_prevalance_summary[2,2]),
                                             sprintf("%0.3f",anxiety_prevalance_summary[3,1]),sprintf("%0.3f",anxiety_prevalance_summary[3,2]),
                                             sprintf("%0.3f",anxiety_prevalance_summary[4,1]),sprintf("%0.3f",anxiety_prevalance_summary[4,2])))

#theme函数设置
theme_bar <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.text = element_text(face = "bold",size = 12),#坐标轴刻度标签加粗
          # axis.ticks = element_line(color='black'),#坐标轴刻度线
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),#去除图例标题
          #legend.justification=c(1,0),#图例在画布的位置(绘图区域外)
          # legend.position=c(0.28, 0.9),#图例在绘图区域的位置
          # legend.position='top',#图例放在顶部
          # legend.direction = "horizontal",#设置图例水平放置
          # legend.spacing.x = unit(2, 'cm'),
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=20)),
          # legend.box = "vertical"
          # legend.box.margin = margin(50, 50, 50, 50, unit = "points")
          # legend.background = element_rect( linetype="solid",colour ="black")
          # legend.margin=margin(0,0,-7,0)#图例与绘图区域边缘的距离
          # legend.box.margin =margin(-10,0,0,0)
    )
  
}

q1<-ggplot(data=plot_data, mapping=aes(x = traits, y = anxiety_prevalence, fill = status))+
  labs(y="Prevalence of anxiety disorder")+    #可以用x,title修改其他名称
  geom_bar(stat="identity",position=position_dodge(0.7),width=0.5)+ #调整柱子宽度
  coord_cartesian(ylim=c(0.00,0.15))+ #设置坐标轴值域
  scale_y_continuous(expand = c(0, 0))+ #消除x轴与绘图区的间隙
  scale_fill_manual(values=alpha(c("#FC4E07","#00AFBB"),0.8))+ #颜色的十六进制代码，或直接用red、blue、green等也可,后面的数字是透明度
  theme_bar()+
  geom_signif(y_position=c(0.065), xmin=c(0.8), xmax=c(1.15),  # 设置显著性说明，y_position是误差线所在y轴位置，xmin和xmax是误差线在x轴位置，可传入多个值
              annotation=c("***"), tip_length=0.02, size=1, textsize = 7,  # 显著性标识；显著性括号下延长度；大小设置；字体大小
              vjust = 0.3)+  # 调整显著性标识和显著性括号之间的距离
  geom_signif(y_position=c(0.138), xmin=c(1.8), xmax=c(2.15),  # 设置显著性说明，y_position是误差线所在y轴位置，xmin和xmax是误差线在x轴位置，可传入多个值
              annotation=c("***"), tip_length=0.02, size=1, textsize = 7,  # 显著性标识；显著性括号下延长度；大小设置；字体大小
              vjust = 0.3)+  # 调整显著性标识和显著性括号之间的距离
  geom_signif(y_position=c(0.140), xmin=c(2.8), xmax=c(3.15),  # 设置显著性说明，y_position是误差线所在y轴位置，xmin和xmax是误差线在x轴位置，可传入多个值
              annotation=c("***"), tip_length=0.02, size=1, textsize = 7,  # 显著性标识；显著性括号下延长度；大小设置；字体大小
              vjust = 0.3)+  # 调整显著性标识和显著性括号之间的距离
  geom_signif(y_position=c(0.117), xmin=c(3.8), xmax=c(4.15),  # 设置显著性说明，y_position是误差线所在y轴位置，xmin和xmax是误差线在x轴位置，可传入多个值
              annotation=c("***"), tip_length=0.02, size=1, textsize = 7,  # 显著性标识；显著性括号下延长度；大小设置；字体大小
              vjust = 0.3)  # 调整显著性标识和显著性括号之间的距离

ggsave('E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/00.UKB_allergy_anxiety_epidemiology/001.anxiety_allergy/column_chi-square.png', plot = q1)
######allergy prevalence in anxiety case and control######
allergy_prevalence_calculation <- function(anxiety_pheno,allergy_pheno)
{
  allergy_prevalance <- data.frame()
  anxiety_temp <- fread(paste0("E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/00.UKB_allergy_anxiety_epidemiology/traits_epidemiology_UKB/",anxiety_pheno,".txt"))
  allergy_temp <- fread(paste0("E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/00.UKB_allergy_anxiety_epidemiology/traits_epidemiology_UKB/",allergy_pheno,".txt"))
  allergy_prevalance[1,1] <- length(intersect(anxiety_temp$eid[which(allergy_temp$allergy_status == 1)],allergy_temp$eid[which(anxiety_temp$anxiety_status == 1)]))/length(which(anxiety_temp$anxiety_status == 1))
  allergy_prevalance[1,2] <- length(intersect(anxiety_temp$eid[which(allergy_temp$allergy_status == 1)],allergy_temp$eid[which(anxiety_temp$anxiety_status == 0)]))/length(which(anxiety_temp$anxiety_status == 0))
  colnames(allergy_prevalance) <- c("anxiety_case_prevalence","anxiety_control_prevalence")
  return(allergy_prevalance)
}

traits <- c("AR_doctor","asthma_doctor","dermatitis_EUR","broad_allergic_disease")
allergy_prevalance_summary <- data.frame()
for(i in 1:length(traits))
{
  if(nrow(allergy_prevalance_summary) == 0)
  {
    allergy_prevalance_summary <- allergy_prevalence_calculation("Anxiety_nerves_generalizedanxietydisorder_professional",traits[i])
  }
  else
  {
    allergy_prevalance_summary <- rbind(allergy_prevalance_summary,allergy_prevalence_calculation("Anxiety_nerves_generalizedanxietydisorder_professional",traits[i]))
  }
}

plot_data <- cbind(data.frame(traits),allergy_prevalance_summary)
plot_data <- data.frame(traits=c(rep("AR",2),rep("Asthma",2),rep("ADE",2),rep("BAD",2)),
                        status=c("Anxiety-Case","Anxiety-Control","Anxiety-Case","Anxiety-Control","Anxiety-Case","Anxiety-Control","Anxiety-Case","Anxiety-Control"),
                        allergy_prevalence=c(allergy_prevalance_summary[1,1],allergy_prevalance_summary[1,2],
                                             allergy_prevalance_summary[2,1],allergy_prevalance_summary[2,2],
                                             allergy_prevalance_summary[3,1],allergy_prevalance_summary[3,2],
                                             allergy_prevalance_summary[4,1],allergy_prevalance_summary[4,2]))




#theme函数设置
theme_bar <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.text = element_text(face = "bold",size = 12),#坐标轴刻度标签加粗
          # axis.ticks = element_line(color='black'),#坐标轴刻度线
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),#去除图例标题
          #legend.justification=c(1,0),#图例在画布的位置(绘图区域外)
          # legend.position=c(0.28, 0.9),#图例在绘图区域的位置
          # legend.position='top',#图例放在顶部
          # legend.direction = "horizontal",#设置图例水平放置
          # legend.spacing.x = unit(2, 'cm'),
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=20)),
          #legend.background = element_rect( linetype="solid",colour ="black")
          # legend.margin=margin(0,0,-7,0)#图例与绘图区域边缘的距离
          # legend.box.margin =margin(-10,0,0,0)
    )
  
}

q1<-ggplot(data=plot_data, mapping=aes(x = traits, y = allergy_prevalence, fill = status))+
  labs(y="Prevalence of allergic disease")+    #可以用x,title修改其他名称
  geom_bar(stat="identity",position=position_dodge(0.7),width=0.5)+ #调整柱子宽度
  coord_cartesian(ylim=c(0.00,0.28))+ #设置坐标轴值域
  scale_y_continuous(expand = c(0, 0))+ #消除x轴与绘图区的间隙
  scale_fill_manual(values=alpha(c("#FC4E07","#00AFBB"),0.8))+ #颜色的十六进制代码，或直接用red、blue、green等也可,后面的数字是透明度
  theme_bar()+
  geom_signif(y_position=c(0.035), xmin=c(0.8), xmax=c(1.15),  # 设置显著性说明，y_position是误差线所在y轴位置，xmin和xmax是误差线在x轴位置，可传入多个值
            annotation=c("***"), tip_length=0.02, size=1, textsize = 7,  # 显著性标识；显著性括号下延长度；大小设置；字体大小
            vjust = 0.3)+  # 调整显著性标识和显著性括号之间的距离
  geom_signif(y_position=c(0.025), xmin=c(1.8), xmax=c(2.15),  # 设置显著性说明，y_position是误差线所在y轴位置，xmin和xmax是误差线在x轴位置，可传入多个值
            annotation=c("***"), tip_length=0.02, size=1, textsize = 7,  # 显著性标识；显著性括号下延长度；大小设置；字体大小
            vjust = 0.3)+  # 调整显著性标识和显著性括号之间的距离
  geom_signif(y_position=c(0.020), xmin=c(2.8), xmax=c(3.15),  # 设置显著性说明，y_position是误差线所在y轴位置，xmin和xmax是误差线在x轴位置，可传入多个值
            annotation=c("***"), tip_length=0.02, size=1, textsize = 7,  # 显著性标识；显著性括号下延长度；大小设置；字体大小
            vjust = 0.3)+  # 调整显著性标识和显著性括号之间的距离
  geom_signif(y_position=c(0.260), xmin=c(3.8), xmax=c(4.15),  # 设置显著性说明，y_position是误差线所在y轴位置，xmin和xmax是误差线在x轴位置，可传入多个值
            annotation=c("***"), tip_length=0.02, size=1, textsize = 7,  # 显著性标识；显著性括号下延长度；大小设置；字体大小
            vjust = 0.3)  # 调整显著性标识和显著性括号之间的距离

ggsave('E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/00.UKB_allergy_anxiety_epidemiology/002.ADE/column_chi-square.png', plot = q1)


#################################logistic visualization: forest plot#####################################
library(forestplot)
library(forestploter)
library(data.table)
library(dplyr)
rm(list = ls())

setwd("E:/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/00.UKB_allergy_anxiety_epidemiology/001.anxiety_allergy/")
logistic_allergy <- fread("logistic_responAnxiety_summary_EUR20240820.txt") %>% as.data.frame()

forest_data <- logistic_allergy[,c(2:5,7)]
#colnames(forest_data)[1:4] <- c("Allergic diseases","OR","95% LCI","95% UCI")
colnames(forest_data)[1:4] <- c("过敏性疾病","OR","95% LCI","95% UCI")
#forest_data$`Allergic diseases` <- c("ADE","AR","Asthma","BAD")
forest_data$过敏性疾病 <- c("特应性皮炎","过敏性鼻炎","哮喘","广义过敏性疾病")
temp <- data.frame()
for(i in 1:nrow(forest_data))
{
  temp[i,1] <- paste0(sprintf("%0.2f",forest_data$OR[i])," (",sprintf("%0.2f",forest_data$`95% LCI`[i]),", ",sprintf("%0.2f",forest_data$`95% UCI`[i]),")")
}
colnames(temp) <- "OR (95% CI)"
forest_data <- cbind(forest_data,temp)
forest_data$`Pr(>|z|)` <- sprintf("%.2e",forest_data$`Pr(>|z|)`)
forest_data$`Pr(>|z|)`[4] <- 0
# 计算标准误差（SE），它在绘图的时候会表示正方形的大小
forest_data$se <- (log(forest_data$`95% UCI`) - log(forest_data$OR))/1.96
forest_data$` ` <- paste(rep("       ", nrow(forest_data)), collapse = " ")

tm <- forest_theme(base_size = 10,           # 设置基础字体大小
                   refline_col = "red4",     # 设置参考线颜色为红色
                   arrow_type = "closed",    # 设置箭头类型为闭合箭头
                   )   

# tm <- forest_theme(
#   base_size = 10,        # 设置文本的基础大小
#   
#   # 设置可信区间的外观
#   ci_pch = 15,           # 可信区间点的形状
#   ci_col = "blue4",    # 可信区间的边框颜色
#   ci_fill = "blue4",      # 可信区间的填充颜色
#   ci_alpha = 0.8,        # 可信区间的透明度
#   ci_lty = 1,            # 可信区间的线型
#   ci_lwd = 1.5,          # 可信区间的线宽
#   ci_Theight = 0.2,      # 设置T字在可信区间末端的高度，默认是NULL
#   
#   # 设置参考线的外观
#   refline_lwd = 0.5,         # 参考线的线宽
#   refline_lty = "dashed",  # 参考线的线型
#   refline_col = "grey20",  # 参考线的颜色
#   
#   # 设置垂直线的外观
#   vertline_lwd = 1,         # 垂直线的线宽，可以添加一条额外的垂直线，如果没有就不显示
#   vertline_lty = "dashed",  # 垂直线的线型
#   vertline_col = "grey20",  # 垂直线的颜色
#   
#   # 设置脚注的字体大小、字体样式和颜色
#   footnote_cex = 0.6,            # 脚注字体大小
#   footnote_fontface = "italic",  # 脚注字体样式
#   footnote_col = "red4"          # 脚注文本的颜色
#   
# )

# 绘制森林图
png(filename = "forest20241027.png",width = 3000,height = 1479,res = 600)
forest(forest_data[,c(1, 8, 6, 5)],   # 选择要在森林图中使用的数据列，这里包括变量名列、患者数量列、绘图要用的空白列和HR（95%CI）列
            est = forest_data$OR,          # 效应值，也就是HR列
            lower = forest_data$`95% LCI`,  # 置信区间下限
            upper = forest_data$`95% UCI`,  # 置信区间上限
            sizes = 0.3,        # 黑框框的大小
            ci_column = 2,             # 在第3列（可信区间列）绘制森林图
            ref_line = 1,              # 添加参考线
            #arrow_lab = c("Protective Factor", "Risk Factor"),  # 箭头标签，用来表示效应方向，如何设置取决于你的样本情况
            xlim = c(0.8, 4.0),          # 设置x轴范围
            ticks_at = c(1),  # 在指定位置添加刻度
            theme = tm)                # 添加自定义主题
            #footnote = "This is the demo data. Please feel free to change\nanything you want.")  # 添加脚注信息
dev.off()
######################info supplement####################
data_dir <- "E:/move/disk_document/Projects_documents/05.allergic_diseases_mental_disorders/01.allergic_diseases_anxiety/00.UKB_allergy_anxiety_epidemiology/traits_epidemiology_UKB/"

trait <- c("Anxiety_nerves_generalizedanxietydisorder_professional_EUR","AR_doctor_EUR","asthma_doctor_EUR","dermatitis_EUR","BAD_EUR")

info <- data.frame()
for(i in 1:length(trait))
{
  data_temp <- fread(paste0(data_dir,trait[i],".txt"))
  colnames(data_temp)[4] <- "allergy_status"
  info[i,1] <- nrow(data_temp)
  info[i,2] <- paste0(nrow(subset(data_temp,allergy_status==1)),"(",sprintf("%0.3f",nrow(subset(data_temp,allergy_status==1))/info[i,1]),")")
  info[i,3] <- paste0(nrow(subset(data_temp,allergy_status==0)),"(",sprintf("%0.3f",nrow(subset(data_temp,allergy_status==0))/info[i,1]),")")
  info[i,4] <- nrow(subset(data_temp,sex==1))
  info[i,5] <- nrow(subset(data_temp,sex==0))
  info[i,6] <- nrow(subset(data_temp,allergy_status==1 & sex==1))
  info[i,7] <- nrow(subset(data_temp,allergy_status==1 & sex==0))
  info[i,8] <- nrow(subset(data_temp,allergy_status==0 & sex==1))
  info[i,9] <- nrow(subset(data_temp,allergy_status==0 & sex==0))
  chi_square_table <- matrix(c(info[i,6],info[i,7],info[i,8],info[i,9]),nrow = 2,ncol = 2)
  info[i,10] <- unname(chisq.test(chi_square_table,correct = F)$p.value)
  info[i,11] <- sprintf("%0.2f",mean(data_temp$age))
  info[i,12] <- sprintf("%0.2f",sd(data_temp$age))
  info[i,13] <- sprintf("%0.2f",mean(subset(data_temp,allergy_status==1)$age))
  info[i,14] <- sprintf("%0.2f",sd(subset(data_temp,allergy_status==1)$age))
  info[i,15] <- sprintf("%0.2f",mean(subset(data_temp,allergy_status==0)$age))
  info[i,16] <- sprintf("%0.2f",sd(subset(data_temp,allergy_status==0)$age))
  info[i,17] <- t.test(subset(data_temp,allergy_status==1)$age,subset(data_temp,allergy_status==0)$age)$p.value
}
fwrite(info,"info_20240820.csv")

anxiety_trait <- "Anxiety_nerves_generalizedanxietydisorder_professional_EUR"
allergy_trait <- c("AR_doctor_EUR","asthma_doctor_EUR","dermatitis_EUR","BAD_EUR")
info_overall_prevalence <- data.frame()
for(i in 1:length(allergy_trait))
{
  anxiety_tmp <- fread(paste0(data_dir,anxiety_trait,".txt"))
  allergy_tmp <- fread(paste0(data_dir,allergy_trait[i],".txt"))
  info_overall_prevalence[1,i] <- sprintf("%0.3f",length(intersect(subset(anxiety_tmp,anxiety_status==1)$eid,allergy_tmp$eid))/nrow(allergy_tmp))
}
