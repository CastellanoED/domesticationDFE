#!/bin/bash
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name=slim_2_vcf
#SBATCH --ntasks=1
#SBATCH --time=02:01:00 
#SBATCH --mem-per-cpu=4gb 
#SBATCH --array=1601-1800
#SBATCH --output=./parallel_jobs/slim_2_vcf-%A_%a.out
#SBATCH --nodes=1

now=$(date +"%T"-"%d"-"%m"-"%Y")
echo "Starting time : $now"


### Break job name into variable params ###
jobsfile=../results/slim_combinations_for_job_array.txt
jobs=$(awk "NR==$SLURM_ARRAY_TASK_ID" $jobsfile)
echo ${jobs}

MIGRATION=$(echo ${jobs} | sed 's/-/ /g' | awk '{print $1}')
echo "Migration:${MIGRATION}"
POSSEL=$(echo ${jobs} | sed 's/-/ /g' | awk '{print $2}')
echo "Positive selection:${POSSEL}"
CHANGE=$(echo ${jobs} | sed 's/-/ /g' | awk '{print $3}')
echo "Change:${CHANGE}"
ITER=$(echo ${jobs} | sed 's/-/ /g' | awk '{print $4}')
echo "Iteration:${ITER}"

### Loading stuff ###
source ~/micromamba/etc/profile.d/mamba.sh 
micromamba activate /home/u12/dcastellano/micromamba/envs/starters/

### Generate VCF-like files from SliM outputs (haplotype information is absent in those VCFs) ###
echo "${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}"
sh extract_mutations_slim.sh ../results/Simulations/sim_${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}_slim_output_sample_file.txt
Rscript --vanilla generate_VCFs_from_Slim_outputs.R ${MIGRATION} ${POSSEL} ${CHANGE} ${ITER} 

now=$(date +"%T"-"%d"-"%m"-"%Y")
echo "Finishing time : $now"


