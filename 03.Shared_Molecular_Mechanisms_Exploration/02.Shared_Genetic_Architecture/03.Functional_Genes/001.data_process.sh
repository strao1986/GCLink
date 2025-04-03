#!/bin/bash
#PBS -N SMR
#PBS -l nodes=cu02:ppn=8
#PBS -e SMR.err
#PBS -o SMR.log
#PBS -q batch
#PBS -V

#trait_anxiety=("anxiety2021" "anxiety2018" "anxiety2017")
trait_allergy=("BAD")
#trait_asthma2018="Asthma2018"

#anxiety_data_dir="/public/jiangjw/02.allergy_anxiety/data/Anxiety"
allergy_data_dir="/public/jiangjw/02.anxiety_ADs/data/Allergic_diseases"

result_dir='/public/jiangjw/02.anxiety_ADs/09.SMR/GWAS_summary_SMR'

#for i in ${trait_anxiety[@]};
#do
#	cat ${anxiety_data_dir}/${i}.clean.txt | cut -f 3-10 | awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$7,$4,$5,$6,$8}' | sed '1d' | sed "1i SNP\tA1\tA2\tfreq\tb\tse\tp\tn" > ${result_dir}/${i}_SMR.clean.txt
#done

for j in ${trait_allergy[@]};
do
	cat ${allergy_data_dir}/${j}.clean.txt | cut -f 3-10 | awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$7,$4,$5,$6,$8}' | sed '1d' | sed "1i SNP\tA1\tA2\tfreq\tb\tse\tp\tn" > ${result_dir}/${j}_SMR.clean.txt
done



#cat ${allergy_data_dir}/Asthma2018.clean.txt | cut -f 3-9  | sed '1d' | sed "1i SNP\tA1\tA2\tb\tse\tp\tn" > ${result_dir}/Asthma2018_SMR.clean.txt



