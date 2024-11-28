#!/bin/bash
#SBATCH --account=
#SBATCH --qos=
#SBATCH --partition=
#SBATCH --job-name="vcf_2_fs"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=07:01:00 
#SBATCH --mem-per-cpu=4gb


for MIGRATION in 0 0.01 
  do      
  for POSSEL in 0 2 20 200
    do
    for CHANGE in 0 0.05 0.25
      do

        ### concatenate here all the VCFs for a given parameter combination
        for ITER in {0..99}
          do
            cat ../results/VCF/sim_${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}_syn.vcf  >> ../results/VCF/sim_${MIGRATION}-${POSSEL}-${CHANGE}_syn.vcf
            cat ../results/VCF/sim_${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}_nsyn.vcf >> ../results/VCF/sim_${MIGRATION}-${POSSEL}-${CHANGE}_nsyn.vcf
            # rm ../results/VCF/sim_${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}_syn.vcf ../results/VCF/sim_${MIGRATION}-${POSSEL}-${CHANGE}-${ITER}_nsyn.vcf
          done

        ### remove all headers except the first one
    		sed '2,${/^#CHROM/d;}' ../results/VCF/sim_${MIGRATION}-${POSSEL}-${CHANGE}_syn.vcf  > ../results/VCF/sim_${MIGRATION}-${POSSEL}-${CHANGE}_syn.one_head.vcf
    		sed '2,${/^#CHROM/d;}' ../results/VCF/sim_${MIGRATION}-${POSSEL}-${CHANGE}_nsyn.vcf > ../results/VCF/sim_${MIGRATION}-${POSSEL}-${CHANGE}_nsyn.one_head.vcf
        rm ../results/VCF/sim_${MIGRATION}-${POSSEL}-${CHANGE}_syn.vcf ../results/VCF/sim_${MIGRATION}-${POSSEL}-${CHANGE}_nsyn.vcf

        ### compute 1d and 2d fs for dadi     
        python dadi_extract_1d_and_2d_sfs.py ${MIGRATION} ${POSSEL} ${CHANGE} syn   #This script computes the 1D-SFS and 2D-SFS, bootstrap in chunks and generates 1D-SFS (for polyDFE) and 2D-SFS (for dadi) using each single ITER as an independent chunk.
        python dadi_extract_1d_and_2d_sfs.py ${MIGRATION} ${POSSEL} ${CHANGE} nsyn 

        ### compute 1d fs for polyDFE
        for BOOT_REPLICATE in {0..99}
          do
            echo "${MIGRATION}-${POSSEL}-${CHANGE}-${BOOT_REPLICATE}"
            Rscript --vanilla adapt_dadi_1D_sfs_to_polyDFE_format.R ${MIGRATION} ${POSSEL} ${CHANGE} ${BOOT_REPLICATE} #This script adapts the bootstrap 1D-SFS output from the previous script to polyDFE input format           
          done

      done
    done
  done


exit 0
