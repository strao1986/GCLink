#!/bin/bash

################################ !! CHANGE !! ##################################
### nomalization GWAS data(.txt) & ref file (.bed)
data_dir='/mnt/e/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/03.GWAS-PW'
ref_dir='/mnt/e/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/reference/nygcresearch-ldetect-data-ac125e47bf7f/EUR/fourier_ls-all.bed' # ATTENTION! ASN/EUR

result_dir='/mnt/e/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/03.GWAS-PW'
################################################################################

#cd ${data_dir}
#list=("anxiety2018_ADE2015.pw.txt" "anxiety2018_ADE2021.pw.txt" "anxiety2018_Asthma2018.pw.txt" "anxiety2018_Asthma2020.pw.txt" "anxiety2018_DiagnosedAR.pw.txt" "anxiety2018_SelfReportAR.pw.txt" "anxiety2021_ADE2015.pw.txt" "anxiety2021_ADE2021.pw.txt" "anxiety2021_Asthma2018.pw.txt" "anxiety2021_Asthma2020.pw.txt" "anxiety2021_DiagnosedAR.pw.txt" "anxiety2021_SelfReportAR.pw.txt") # In cycling, use `ls`; for single run, use manually assignment --23.12.27

list=("anxiety2018_BAD.pw.txt" "anxiety2021_BAD.pw.txt")
### where GNU & sBOOST located (dependicies of GWAS-PW)
## The program should run in this path where GNU & sBOOST located (dependicies of GWAS-PW)
##lib_dir='/public/chenxt/miniconda3/envs/gwas-pw/lib'
##cd ${lib_dir}

for i in ${list[@]};do 
    # Filter ".dist" files.
	filename=`echo ${i%%.pw.*}` #extracts the base filename by removing the ".pw.*" part from the filename. The double percentage sign (%%) is used for greedy matching from the end of the string.
    ## According to the name of input file, the phenotype before "_" is Phenotype1.
    mkdir -p ${result_dir}/${filename}
	pheno1=`echo ${filename%%_*}`
	pheno2=`echo ${filename##*_}`
    ./gwas-pw-0.21/src/gwas-pw \
    -i ${data_dir}/${i} \
	-bed ${ref_dir} \
	-phenos ${pheno1} ${pheno2} \
	-o ${result_dir}/${filename}/${pheno1}_${pheno2}
	#-o ${result_dir}/${filename}/${pheno1}_${pheno2} >  ${result_dir}/${filename}/gwas-pw_${pheno1}_${pheno2}.log 2>&1  
# 	echo "nohup gwas-pw \
# 		-i ${data_dir}/${i} \
# 		-bed ${ref_dir} \
# 		-phenos ${pheno1} ${pheno2} \
# 		-o ${result_dir}/${filename}/${pheno1}_${pheno2} >  ${result_dir}/${filename}/gwas-pw_${pheno1}_${pheno2}.log 2>&1  & " | bash
done
