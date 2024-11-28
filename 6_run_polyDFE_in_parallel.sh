#!/bin/bash
#SBATCH --account=
#SBATCH --qos=
#SBATCH --partition=
#SBATCH --job-name=polyDFE_domestication
#SBATCH --ntasks=1
#SBATCH --time=96:00:00 
#SBATCH --mem-per-cpu=4gb 
#SBATCH --array=1-3000
#SBATCH --output=./parallel_jobs/polyDFE_domestication-%A_%a.out
#SBATCH --nodes=1

now=$(date +"%T"-"%d"-"%m"-"%Y")
echo "Starting time : $now"

### Break job name into variable params ###
jobsfile=../results/polyDFE_jobs_list.txt 
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
MODEL=$(echo ${jobs} | sed 's/-/ /g' | awk '{print $5}')
echo "polyDFE Model:${MODEL}"


./polyDFE-2.0-linux-64-bit -d ../results/SFS/for_polyDFE/sim_${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}-wild.fs:../results/SFS/for_polyDFE/sim_${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}-domesticated.fs -m C -t -o bfgs -i init_model_BandC.txt ${MODEL}  -r range_model_BandC.txt 0 -v 50 -b params.basinhop 1 -w -l -1 > ../results/polyDFE_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}_${MODEL}_polyDFE.txt  

now=$(date +"%T"-"%d"-"%m"-"%Y")
echo "Finishing time : $now"

exit 0

