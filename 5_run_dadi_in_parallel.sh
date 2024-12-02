#!/bin/bash
#SBATCH --account=
#SBATCH --qos=
#SBATCH --partition=
#SBATCH --job-name=dadi_domestication
#SBATCH --ntasks=1 #94 when generating the caches, otherwise 1 task
#SBATCH --time=04:00:00 #3 days when generating caches, otherwise 4 hours
#SBATCH --mem-per-cpu=4gb 
#SBATCH --array=1-24 ##first run from 1 to 24, to infer demographic parameters and cache for the first time with the original data, then run from 25 to 2400
#SBATCH --output=./parallel_jobs/dadi_domestication-%A_%a.out
#SBATCH --nodes=1


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


if [ ${ITER} == 0 ]
then
	### The demography is inferred using the original, non-bootstrapped, 2D-SFS. 
	echo "Inferring domestication demography and generating cache using original data"
	python dadi_infer_flexible_demography_2d.py ${MIGRATION} ${POSSEL} ${CHANGE}
	wait

	Rscript best_fit.R ${MIGRATION} ${POSSEL} ${CHANGE} flexible

	python dadi_plot_demographic_fit.py   ${MIGRATION} ${POSSEL} ${CHANGE} complex 
	python dadi_uncertainty_demography.py ${MIGRATION} ${POSSEL} ${CHANGE} # Here the bootstrapped synonymous 2D-SFS are used for the Godambe uncertainty inference.
	python dadi_generate_cache.py         ${MIGRATION} ${POSSEL} ${CHANGE}

	echo "Inferring domestication 2D-DFE using original data"
	### Compute loglk of each scenario under 4 positive DFEs Sb: 0, 1, 10 and 100.
	### Compute loglk of each scenario where ppos_wild == pchange_pos (true model) vs ppos_wild != pchange_pos (current agnostic model) 
	python dadi_test_positive_DFE.py ${MIGRATION} ${POSSEL} ${CHANGE}

	sort -k1,1n ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_fits_with_Sb=0_two_pb.txt | tail -n 1 > ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_with_Sb=0_best_fit.txt
	sort -k1,1n ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_fits_with_Sb=2_two_pb.txt | tail -n 1 > ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_with_Sb=2_best_fit.txt
	sort -k1,1n ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_fits_with_Sb=20_two_pb.txt | tail -n 1 > ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_with_Sb=20_best_fit.txt
	sort -k1,1n ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_fits_with_Sb=200_two_pb.txt | tail -n 1 > ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_with_Sb=200_best_fit.txt
	                      
	python dadi_plot_dfe_fit.py ${MIGRATION} ${POSSEL} ${CHANGE} 0 &
	python dadi_plot_dfe_fit.py ${MIGRATION} ${POSSEL} ${CHANGE} 2 &
	python dadi_plot_dfe_fit.py ${MIGRATION} ${POSSEL} ${CHANGE} 20 &
	python dadi_plot_dfe_fit.py ${MIGRATION} ${POSSEL} ${CHANGE} 200 &
	wait
fi

echo "Inferring domestication 2D-DFE across bootstrap replicates"
python dadi_infer_dfe.py  ${MIGRATION} ${POSSEL} ${CHANGE} ${ITER} # The uncertainty of the DFE is inferred using the tradional bootstrap "in chunks" approach. 

now=$(date +"%T"-"%d"-"%m"-"%Y")
echo "Finishing time : $now"

exit 0

