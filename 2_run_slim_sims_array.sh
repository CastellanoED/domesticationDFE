#!/bin/bash
#SBATCH --account=
#SBATCH --qos=
#SBATCH --partition=
#SBATCH --job-name=domestication_sims
#SBATCH --ntasks=1
#SBATCH --time=120:01:00 
#SBATCH --mem-per-cpu=4gb 
#SBATCH --array=1-2400
#SBATCH --output=./parallel_jobs/domestication_sims-%A_%a.out
#SBATCH --nodes=1

now=$(date +"%T"-"%d"-"%m"-"%Y")
echo "Starting time : $now"

### Fix params ###
NLOCI=10000
SIZEG=120
NEPOP=5000
NEBOT=25
MEAN_DEL=-200
let TBOT=200*${NEPOP}/10000
let TEND=1800*${NEPOP}/10000
let NE1D=${NEPOP}/${NEBOT}

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

if [ ${POSSEL}==0 ] #this is a sanity check where all nsyn are deleterious
then 
  DEL=1
  MEAN_DEL=-200
fi

if [ ${POSSEL}==2 ]
then 
  DEL=0.9
fi

if [ ${POSSEL}==20 ]
then 
  DEL=0.99
fi

if [ ${POSSEL}==200 ]
then 
  DEL=0.999
fi

echo "${NLOCI}-${SIZEG}-${NEPOP}-${NEBOT}-${TBOT}-${TEND}-${NE1D}-${ITER}"
echo "${MIGRATION}-${POSSEL}-${CHANGE}-${DEL}-${ITER}"

### Run SliM ###
slim -t -m -d "seed=$RANDOM" -d "nloci=$NLOCI" -d "nchrom=1" -d "size_gene=$SIZEG" \
     -d "Ne_pop=$NEPOP" \
     -d "Ne1d=$NE1D" \
     -d "t0=10*Ne_pop" \
     -d "t1=1" \
     -d "tbot=$TBOT" \
     -d "tend=$TEND" \
     -d "migration_rate=$MIGRATION" \
     -d "mutation_rate=1.25e-7*integerDiv(10000,Ne_pop)" \
     -d "maprec='flat'" -d "recrate_within=0.6*mutation_rate" \
	   -d "sample_size1=40" \
     -d "s_mean_beneficial=$POSSEL/(2*Ne_pop)" \
     -d "s_mean_deleterious=$MEAN_DEL/(2*Ne_pop)" -d "shape_deleterious=0.3" -d "h=0.5" -d "prop_del_anc=$DEL" -d "change_prop=$CHANGE" -d "prop_del_new=$DEL" \
     -d "Rfile_ref='./do_REF_A.R'" \
     -d "fileREF_fa='../results/fastas/sim_${jobs}_REF_A.fa'" \
     -d "Rfpm_file='./calculate_fitness_position_matrix.R'" \
     -d "file_fm='../results/fitness_position_matrix/sim_${jobs}_fitness_position_matrix.txt'" \
     -d "file_output1='../results/Simulations/sim_${jobs}_slim_output_sample_file.txt'" \
     -d "file_log1='../results/Simulations/sim_${jobs}_fitness_out_wild_before_dom.txt'" \
     -d "file_log2='../results/Simulations/sim_${jobs}_fitness_out_wild_dom_end.txt'" \
     ./slim_code_mod4_NEW.slim > ../results/Simulations/sim_${jobs}_exit.txt 2> ../results/Simulations/sim_${jobs}_time.txt 

now=$(date +"%T"-"%d"-"%m"-"%Y")
echo "Finishing time : $now"


