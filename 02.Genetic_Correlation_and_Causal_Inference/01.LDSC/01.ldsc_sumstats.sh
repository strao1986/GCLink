#/bin/bash
#PBS -S /bin/bash
#PBS -N ldsc
#PBS -l nodes=cu04:ppn=5
#PBS -e ldsc.err
#PBS -o ldsc.log
#PBS -q batch
#PBS -V

export OMP_NUM_THREADS=5
export OPENBLAS_NUM_THREADS=5
export MKL_NUM_THREADS=5
export VECLIB_MAXIMUM_THREADS=5
export NUMEXPR_NUM_THREADS=5

source ~/miniconda3/bin/activate ldsc

########################## CHANGE ##########################
ref_dir='/public/jiangjw/02.anxiety_ADs/reference/eur_w_ld_chr'
#data_path='/public/jiangjw/02.anxiety_ADs/data/Allergic_diseases'
TN_path='/public/jiangjw/05.UkbGWAS/trigeminal_neuralgia'
MDD_path='/public/jiangjw/06.MDD_ADs/data/EUR'

traits_TN=("TrigeminalNeuralgia")
traits_MDD=("MDD")

result_dir='/public/jiangjw/05.UkbGWAS/trigeminal_neuralgia/TN_MDD/LDSC'

run_dir='/public/jiangjw/02.anxiety_ADs/02.LDSC/ldsc-master'
############################################################

## 1. Munge Data
# Length of array: # before arry, and [@] after arry(#array[@])
## --N [number] set the sample size of data.
## --N-col [colnamne of sample size ] specify the colname of Sample size.
## --a1 [A1] Name of A1 column
for i in ${traits_TN[@]}
do 
	python ${run_dir}/munge_sumstats.py \
	--sumstats ${TN_path}/${i}.clean.txt \
	--N-col N \
	--out ${result_dir}/${i} \
	--merge-alleles ${ref_dir}/w_hm3.snplist \
	--chunksize 500000 # if this option not set, the program may not run.
done

for j in ${traits_MDD[@]}
do
        python ${run_dir}/munge_sumstats.py \
        --sumstats ${MDD_path}/${j}.clean.txt \
        --N-col N \
        --out ${result_dir}/${j} \
        --merge-alleles ${ref_dir}/w_hm3.snplist \
        --chunksize 500000 # if this option not set, the program may not run.
done

