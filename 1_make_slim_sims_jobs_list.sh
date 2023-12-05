#!/bin/bash
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="Slim_params"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:01:00 
#SBATCH --mem-per-cpu=1gb

rm ../results/slim_combinations_for_job_array.txt
rm ../results/dadi_combinations_for_job_array.txt

for ITER in {1..100}
  do
    for MIGRATION in 0 0.01
      do      
        for POSSEL in 2 20 200
          do
            for CHANGE in 0 0.05 0.25
              do
                echo "${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}" >> ../results/slim_combinations_for_job_array.txt
                echo "${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}" >> ../results/dadi_combinations_for_job_array.txt

              done
          done
      done
  done

### Adding some extra neutral simulations as sanity check ###
for ITER in {1..10}
  do
    for MIGRATION in 0 0.01
      do      
        for POSSEL in 0
          do
            for CHANGE in 0 
              do
                echo "${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}" >> ../results/slim_combinations_for_job_array.txt
                echo "${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}" >> ../results/dadi_combinations_for_job_array.txt

              done
          done
      done
  done


exit 0

