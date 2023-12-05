#!/bin/bash
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name=dadi_domestication
#SBATCH --ntasks=1 #94 when generating the caches, otherwise 1 task
#SBATCH --time=24:00:00 #3 days when generating caches, otherwise 4 hours
#SBATCH --mem-per-cpu=4gb 
#SBATCH --array=1783-1800 ##first run from 1 to 18, to infer demographic parameters and cache for the first time with the original data, then run from 19 to 1800
#SBATCH --output=./parallel_jobs/dadi_domestication-%A_%a.out
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


if [ ${ITER} == 1 ]
then
	echo "Inferring domestication demography and generating cache"
	### python dadi_infer_demography_2d.py         ${MIGRATION} ${POSSEL} ${CHANGE} &
	### python dadi_infer_complex_demography_2d.py ${MIGRATION} ${POSSEL} ${CHANGE}
	python dadi_infer_flexible_demography_2d.py ${MIGRATION} ${POSSEL} ${CHANGE}
	wait

	micromamba activate /home/u12/dcastellano/micromamba/envs/starters/
	### Rscript best_fit.R ${MIGRATION} ${POSSEL} ${CHANGE} true
	### Rscript best_fit.R ${MIGRATION} ${POSSEL} ${CHANGE} complex
	Rscript best_fit.R ${MIGRATION} ${POSSEL} ${CHANGE} flexible

	micromamba activate /home/u12/dcastellano/micromamba/envs/dadi/
	python dadi_plot_demographic_fit.py   ${MIGRATION} ${POSSEL} ${CHANGE} complex 
	python dadi_uncertainty_demography.py ${MIGRATION} ${POSSEL} ${CHANGE}
	python dadi_generate_cache.py         ${MIGRATION} ${POSSEL} ${CHANGE}
fi

if [ ${ITER} == 100 ]
then
	ITER=0
fi

echo "Inferring domestication 2D-DFE across bootstrap replicates"
python dadi_infer_dfe.py  ${MIGRATION} ${POSSEL} ${CHANGE} ${ITER}


now=$(date +"%T"-"%d"-"%m"-"%Y")
echo "Finishing time : $now"

exit 0

