#!/bin/bash
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name=dadi_test_positive_DFE
#SBATCH --ntasks=1 
#SBATCH --time=24:00:00 
#SBATCH --mem-per-cpu=4gb 
#SBATCH --array=1-18 
#SBATCH --output=./parallel_jobs/dadi_test_positive_DFE-%A_%a.out
#SBATCH --nodes=1

source ~/micromamba/etc/profile.d/mamba.sh 
micromamba activate /home/u12/dcastellano/micromamba/envs/dadi/

now=$(date +"%T"-"%d"-"%m"-"%Y")
echo "Starting time : $now"

### Break job name into variable params ###
jobsfile=../results/dadi_combinations_for_job_array.txt
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

### Original data analysis
### Compute loglk of each scenario under 4 positive DFEs Sb: 0, 1, 10 and 100.
### Compute loglk of each scenario where ppos_wild == pchange_pos (true model) vs ppos_wild != pchange_pos (current model) 

python dadi_test_positive_DFE.py ${MIGRATION} ${POSSEL} ${CHANGE}

now=$(date +"%T"-"%d"-"%m"-"%Y")
echo "Finishing time : $now"

exit 0

