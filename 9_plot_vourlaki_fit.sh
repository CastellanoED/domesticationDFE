#!/bin/bash
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name=vourlaki_fit_plots
#SBATCH --ntasks=1 
#SBATCH --time=00:30:00 
#SBATCH --mem-per-cpu=4gb 
#SBATCH --array=1-18 
#SBATCH --output=./parallel_jobs/vourlaki_fit_plots-%A_%a.out
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
sort -k1,1n ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_fits_with_Sb=0_two_pb.txt | tail -n 1 > ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_with_Sb=0_best_fit.txt
sort -k1,1n ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_fits_with_Sb=2_two_pb.txt | tail -n 1 > ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_with_Sb=2_best_fit.txt
sort -k1,1n ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_fits_with_Sb=20_two_pb.txt | tail -n 1 > ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_with_Sb=20_best_fit.txt
sort -k1,1n ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_fits_with_Sb=200_two_pb.txt | tail -n 1 > ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_with_Sb=200_best_fit.txt
                      
python dadi_plot_dfe_fit.py ${MIGRATION} ${POSSEL} ${CHANGE} 0
python dadi_plot_dfe_fit.py ${MIGRATION} ${POSSEL} ${CHANGE} 2
python dadi_plot_dfe_fit.py ${MIGRATION} ${POSSEL} ${CHANGE} 20
python dadi_plot_dfe_fit.py ${MIGRATION} ${POSSEL} ${CHANGE} 200

now=$(date +"%T"-"%d"-"%m"-"%Y")
echo "Finishing time : $now"

exit 0

