#!/bin/bash
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="wrap_up"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:01:00 
#SBATCH --mem-per-cpu=1gb


### Wrap up demographic inference ###
rm ../results/col1.txt ../results/col2.txt ../results/col3.txt  ../results/col4.txt

for POSSEL in 2 20 200  
  do
    for MIGRATION in 0 0.01
      do
        for CHANGE in 0 0.05 0.25
          do
              echo "${MIGRATION}-${POSSEL}-${CHANGE}" 
              echo "${MIGRATION}-${POSSEL}-${CHANGE}" >> ../results/col1.txt
              
              echo "${MIGRATION}-${POSSEL}-${CHANGE}" >> ../results/col3.txt
              echo "${MIGRATION}-${POSSEL}-${CHANGE}" >> ../results/col3.txt
              echo "${MIGRATION}-${POSSEL}-${CHANGE}" >> ../results/col3.txt
              echo "${MIGRATION}-${POSSEL}-${CHANGE}" >> ../results/col3.txt
              echo "${MIGRATION}-${POSSEL}-${CHANGE}" >> ../results/col3.txt
              echo "${MIGRATION}-${POSSEL}-${CHANGE}" >> ../results/col3.txt

              sort -k1,1n ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}.2d.flexible_dem_syn_best_fit.txt >> ../results/col2.txt
              cat ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}_demographic_confidence_intervals.txt >> ../results/col4.txt
          done
      done
  done

paste ../results/col1.txt ../results/col2.txt > ../results/data_frames/flexible_dem_best_fits.txt
paste ../results/col3.txt ../results/col4.txt > ../results/data_frames/flexible_dem_Godambe_CI.txt

## Wrap up DFE inference ###
rm ../results/col1.txt ../results/col2.txt
rm ../results/col3.txt ../results/col4.txt

for NEWPOSSEL in 0 2 20 200
  do   
    for pb in two_pb uniq_pb #or uniq_pb for the original_fs
      do
        for POSSEL in 2 20 200  
          do
            for MIGRATION in 0 0.01
              do
                for CHANGE in 0 0.05 0.25
                  do
                      echo "${MIGRATION}-${POSSEL}-${CHANGE}-${NEWPOSSEL}-${pb}" 
                      echo "${MIGRATION}-${POSSEL}-${CHANGE}-${NEWPOSSEL}-${pb}" >> ../results/col1.txt

                      sort -k1,1n ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_fits_with_Sb=${NEWPOSSEL}_${pb}.txt | tail -n 1 >> ../results/col2.txt

                      sort -k1,1n ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_fits_with_Sb=${NEWPOSSEL}_${pb}.txt | tail -n 1 > ../results/dadi_outputs/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_with_Sb=${NEWPOSSEL}_best_fit.txt
                      
                      for ITER in `seq 0 99`
                        do

                          echo "${ITER}-${MIGRATION}-${POSSEL}-${CHANGE}-${NEWPOSSEL}-${pb}" 
                          echo "${ITER}-${MIGRATION}-${POSSEL}-${CHANGE}-${NEWPOSSEL}-${pb}" >> ../results/col3.txt

                          sort -k1,1n ../results/dadi_outputs/Bootstraps/sim_${MIGRATION}-${POSSEL}-${CHANGE}-nsyn_boot.${ITER}_fits_with_Sb=${NEWPOSSEL}_${pb}.txt | tail -n 1 >> ../results/col4.txt
                        done
                  done
              done
          done
      done  
  done

paste ../results/col1.txt ../results/col2.txt > ../results/data_frames/vourlaki_DFE_original_fs_best_fits.txt
paste ../results/col3.txt ../results/col4.txt > ../results/data_frames/vourlaki_DFE_bootstra_fs_best_fits.txt


### Wrap up polyDFE outputs ###

source ~/micromamba/etc/profile.d/mamba.sh 
micromamba activate /home/u12/dcastellano/micromamba/envs/starters/

rm ../results/data_frames/polyDFE_lrt_replicate.txt ../results/data_frames/polyDFE_DFE_replicates_allmodels.txt ../results/data_frames/polyDFE_parameters_replicates_allmodels.txt ../results/data_frames/polyDFE_alpha_replicates_allmodels.txt

for MIGRATION in 0 0.01
  do      
    for POSSEL in 2 20 200
      do
        for CHANGE in 0 0.05 0.25
          do
            for ITER in {1..99} #this should be from 0 to 99
              do
                echo "${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}" 
                Rscript --vanilla polyDFE_output_analyses.R ${MIGRATION} ${POSSEL} ${CHANGE} ${ITER}
              done
          done
      done
  done

Rscript --vanilla polyDFE_output_figures.R



exit 0
