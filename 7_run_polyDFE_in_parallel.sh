#!/bin/bash
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name=polyDFE_domestication
#SBATCH --ntasks=1
#SBATCH --time=96:00:00 
#SBATCH --mem-per-cpu=4gb 
#SBATCH --array=8001-9000
#SBATCH --output=./parallel_jobs/polyDFE_domestication-%A_%a.out
#SBATCH --nodes=1

now=$(date +"%T"-"%d"-"%m"-"%Y")
echo "Starting time : $now"

### Break job name into variable params ###
jobsfile=../results/polyDFE_jobs_list.txt ## the true sfs start by 0 and end at 99, however originally in this list they went from 1 to 100. CORRECTED by replacing 100 by 0.
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

if [ ${ITER} == 100 ]
then
	ITER=0
fi

./polyDFE-2.0-linux-64-bit -d ../results/SFS/for_polyDFE/sim_${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}-wild.fs:../results/SFS/for_polyDFE/sim_${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}-domesticated.fs -m C -t -o bfgs -i ./polyDFEv2.0/custom_input/init_model_BandC.txt ${MODEL}  -r ./polyDFEv2.0/custom_input/range_model_BandC.txt 0 -v 50 -b ./polyDFEv2.0/custom_input/params.basinhop 1 -w -l -1 > ../results/polyDFE_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}_${MODEL}_polyDFE.txt  

now=$(date +"%T"-"%d"-"%m"-"%Y")
echo "Finishing time : $now"

exit 0

# ./polyDFE-2.0-linux-64-bit -d ../results/SFS/for_polyDFE/sim_0-2-0-1-wild.fs:../results/SFS/for_polyDFE/sim_0-2-0-1-domesticated.fs -m C -t -o bfgs -i ./polyDFEv2.0/custom_input/init_model_BandC.txt 20  -r ./polyDFEv2.0/custom_input/range_model_BandC.txt 0 -v 50 -b ./polyDFEv2.0/custom_input/params.basinhop 1 -w -l -1 > ../results/polyDFE_outputs/sim_0-2-0-1_20_polyDFE.txt  
