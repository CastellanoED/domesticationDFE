#!/bin/bash
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="polyDFE_jobs_list"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:01:00 
#SBATCH --mem-per-cpu=1gb

rm ../results/polyDFE_jobs_list.txt

for MIGRATION in 0 0.01
  do      
    for POSSEL in 2 20 200
      do
        for CHANGE in 0 0.05 0.25
          do
            for ITER in {1..100} #this should be from 0 to 99
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
