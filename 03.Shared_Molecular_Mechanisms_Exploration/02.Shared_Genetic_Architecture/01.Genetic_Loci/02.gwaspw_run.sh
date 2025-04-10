#!/bin/bash

GWASPW_input_dir='/mnt/e/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/Github_GCLink/03.Shared_Molecular_Mechanisms_Exploration/02.Shared_Genetic_Architecture/01.Genetic_Loci' #modify into your path
reference_dir='/mnt/e/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/Github_GCLink/03.Shared_Molecular_Mechanisms_Exploration/02.Shared_Genetic_Architecture/01.Genetic_Loci/nygcresearch-ldetect-data-ac125e47bf7f/EUR/fourier_ls-all.bed' #modify into your path
result_dir='/mnt/e/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/Github_GCLink/03.Shared_Molecular_Mechanisms_Exploration/02.Shared_Genetic_Architecture/01.Genetic_Loci' #modify into your path
run_dir='/mnt/e/Projects_documents/05.allergic_diseases_mental_disorders/01.anxiety_ADs/03.GWAS-PW/gwas-pw-0.21/src' #modify into your path

pheno_list=("Anxiety2021_DiagnosedAR.pw.txt")

#### GWAS-PW analysis ####
for i in ${pheno_list[@]};do 
	filename=`echo ${i%%.pw.*}`
	pheno1=`echo ${filename%%_*}`
	pheno2=`echo ${filename##*_}`
	${run_dir}/gwas-pw \
	-i ${GWASPW_input_dir}/${i} \
	-bed ${reference_dir} \
	-phenos ${pheno1} ${pheno2} \
	-o ${result_dir}/${pheno1}_${pheno2}
done
