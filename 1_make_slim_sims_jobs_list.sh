#!/bin/bash
#SBATCH --account=
#SBATCH --qos=
#SBATCH --partition=
#SBATCH --job-name="Slim_params"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:01:00 
#SBATCH --mem-per-cpu=1gb

rm ../results/slim_combinations_for_job_array.txt
rm ../results/dadi_combinations_for_job_array.txt
rm ../results/polyDFE_jobs_list.txt


for ITER in {0..99}
  do
    for MIGRATION in 0 0.01
      do      
        for POSSEL in 0 2 20 200
          do
            for CHANGE in 0 0.05 0.25
              do
                echo "${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}" >> ../results/slim_combinations_for_job_array.txt
                echo "${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}" >> ../results/dadi_combinations_for_job_array.txt

              done
          done
      done
  done

for MIGRATION in 0 0.01
  do      
    for POSSEL in 0 2 20 200
      do
        for CHANGE in 0 0.05 0.25
          do
            for ITER in {0..99} 
              do
                for MODEL in 30 1 2 10 20
                  do
                      echo "${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}-${MODEL}" >> ../results/polyDFE_jobs_list.txt
                  done
              done
          done
      done
  done




exit 0

