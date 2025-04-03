#!/bin/bash
#PBS -S /bin/bash
#PBS -N ldsc
#PBS -l nodes=cu03:ppn=1
#PBS -e ldsc.err
#PBS -o ldsc.log
#PBS -q batch
#PBS -V

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

source ~/miniconda3/bin/activate ldsc

########################## CHANGE ##########################
LDscore_ref_path='/public/jiangjw/reference/eur_w_ld_chr/'
result_dir='/public/jiangjw/02.anxiety_ADs_EUR/02.LDSC'
ldsc_result_dir='/public/jiangjw/02.anxiety_ADs_EUR/02.LDSC'
run_dir='/public/jiangjw/02.anxiety_ADs_EUR/02.LDSC/ldsc-master'
############################################################

### 2. LD Score Regression
## 1) LD Regression
# ldsc.py \
#    --rg mdd_eas.sumstats.gz,sa.sumstats.gz \
#    --ref-ld-chr ${ref_path} \
#    --w-ld-chr ${ref_path} \
#    --out ${ldsc_result_dir}/mddEAS_sa

Anxiety=("Anxiety2021" "Anxiety2018")
#Allergic_disease=("ADE2021" "ADE2015" "AR_Diagnose" "AR_SelfReport" "Asthma2020" "Asthma2018")
#IL=("IL-1a" "IL-1b" "IL-4" "IL-6")
Allergic_disease="BAD"
 
for i in ${Anxiety[@]}
do
	for j in ${Allergic_disease[@]}
	do  
		python ${run_dir}/ldsc.py \
   		--rg ${ldsc_result_dir}/${i}.sumstats.gz,${ldsc_result_dir}/${j}.sumstats.gz \
   		--ref-ld-chr ${LDscore_ref_path} \
   		--w-ld-chr ${LDscore_ref_path} \
   		--out ${result_dir}/ldsc_${i}_${j}
	done
done

#for k in ${Allergic_disease[@]}
#do
#        for l in ${IL[@]}
#        do
#                python ${run_dir}/ldsc.py \
#                --rg ${ldsc_result_dir}/${k}.sumstats.gz,${ldsc_result_dir}/${l}.sumstats.gz \
#                --ref-ld-chr ${ref_path} \
#                --w-ld-chr ${ref_path} \
#                --out ${result_dir}/ldsc_${k}_${l}
#        done
#done
# ldsc.py \
#    --rg mdd.sumstats.gz,sa_fzscore2_z.sumstats.gz \
#    --ref-ld-chr ${ref_dir}/eas_ldscores/ \
#    --w-ld-chr ${ref_dir}/eas_ldscores/ \
#    --out mdd_safzscore2

## 2）Calculate heritability
#for i in ${Allergic_disease[@]}
#do
#	python ${run_dir}/ldsc.py \
#	--h2 ${ldsc_result_dir}/${i}.sumstats.gz \
#	--ref-ld-chr ${ref_path} \
#	--w-ld-chr ${ref_path} \
#	--out ${ldsc_result_dir}/${i}_h2
#done

#for j in ${Anxiety[@]}
#do
#        python ${run_dir}/ldsc.py \
#        --h2 ${ldsc_result_dir}/${j}.sumstats.gz \
#        --ref-ld-chr ${ref_path} \
#        --w-ld-chr ${ref_path} \
#        --out ${ldsc_result_dir}/${j}_h2
#done

# ldsc.py \
#    --h2 sa.sumstats.gz \
#    --ref-ld-chr ${ref_path} \
#    --w-ld-chr ${ref_path} \
#    --out ${ldsc_result_dir}/sa_h2
