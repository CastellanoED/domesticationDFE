## Detection of Domestication Signals through the Analysis of the Full Distribution of Fitness Effects using Simulations

We use forward-in-time simulations, which let us study the domestication process under different demographic and selective models. This helps us understand how well we can detect the selective differences between these two populations (domesticated and wild). We've tried out lots of different combinations of genetic architectures and selective effects (both beneficial and deleterious). Some have a relatively small number of loci that change their selective effects, while others consider polygenic adaptation (where fitness is seen as a trait) where many loci have divergent selective effects. We're excited to announce the release of a new methodology that includes an additional parameter crucial to distinguishing populations in processes of rapid selective change. This parameter allows us to see how the selective effects of a fraction of the existing variants can change (from deleterious to beneficial and vice versa) in the domesticated population. This method is really clever because it can work out all the DFE parameters for both the wild and domesticated populations at the same time. 

More details can be found in:

    https://doi.org/10.1101/2022.08.24.505198

If using this software or simulations in a publication, please cite the above article.

## Contact

David Castellano, dcastellano at arizona . edu, for bug reports and questions about how this pipeline works.
Sebastian Ramos-Onsins, sebastian . ramos at cragenomica . es, for questions about the SLIM simulation recipe (slim_code_mod4_NEW.slim) and how we assign selection coefficients to different sites and populations in the genome (calculate_fitness_position_matrix.R). 
Ryan Gutenkunst, rgutenk at arizona . edu, for questions about the new functions to infer demography, the DFE and the fraction of mutations with divergent selection (domestication_new_dadi_functions.py)

## Installation

Clone the repository:
   ```bash
   git clone https://github.com/CastellanoED/domesticationDFE.git
   ```

## Requirements

- SLiM4 (https://github.com/MesserLab/SLiM)
- dadi (https://bitbucket.org/gutenkunstlab/dadi/)
- polyDFE (https://github.com/paula-tataru/polyDFE)
- R
- git

## Steps and Scripts Overview

Hello there! Before we get started, let's set up some working directories to store all the simulations, DFE inferences, and other goodies. In your preferred local directory, just type:

```bash
mkdir -p domesticationDFE/scripts \
         domesticationDFE/results/dadi_outputs/Bootstraps  \
         domesticationDFE/results/data_frames  \
         domesticationDFE/results/fastas  \
         domesticationDFE/results/fitness_position_matrix  \
         domesticationDFE/results/polyDFE_outputs  \
         domesticationDFE/results/SFS/by_dadi  \
         domesticationDFE/results/SFS/for_polyDFE  \
         domesticationDFE/results/Simulations  \
         domesticationDFE/results/VCF  \
         domesticationDFE/plots/by_dadi \
         domesticationDFE/plots/dadi \
         domesticationDFE/plots/polyDFE
```
If you're using this GitHub repo (which has the latest version of this pipeline), you'll need to move its contents to the domesticationDFE/scripts/ folder and work from there. Let us know if you need any help with this! On the other hand, if you're using the ZENODO stable repo (which contains the version just after peer review), you're all set! No need to do anything else.

The bash scripts are set up in a specific order to make sure that each task is done right: 


```sh 1_make_slim_sims_jobs_list.sh
```

This first bash script simply generates the list of all the 2,400 names/jobs/files we want to simulate and analyse. This is stored in domesticationDFE/results/.


```sh 2_run_slim_sims_array.sh
```

This is where we actually simulate the 24 scenarios (100 replicates for each), all simulation results are stored in domesticationDFE/results/Simulations directory. This is a VERY COMPUTATIONALLY INTENSIVE STEP. Just a quick note to let you know that at ZENODO (), we've made it super easy for you to avoid running all the simulations. 


```sh 3_process_slim_outputs_VCF.sh
```

This third bash script takes the SLiM output from the previous step and converts it into VCF-like files. Note that since we are only interested in allele frequencies, haplotype information is missing in these VCFs, genotypes are not phased, but they could potentially be if the information from the "Genomes:" section in "output_sample_file.txt" is used.


```sh 4_process_slim_outputs_SFS.sh
```

The fourth bash script converts the VCF for each parameter combination into 100 boostrapped SFS for synonymous and nonsynonymous polymorphisms. Each simulation iteration is treated as an independent chromosome chunk. The 1D and 2D-SFS are then converted to polyDFE and dadi format respectively. We've provided the contents of the domesticationDFE/results/SFS/ folder. We hope this makes your life a little easier! If you're a user interested in evaluating your own software in our simulations, we'd love to have you! You can find the 1D-SFS and 2D-SFS in the domesticationDFE/results/SFS/ folder.

# Note to self: Provide ../results/Wild_Domesticated_individuals.txt file

```sh 5_run_dadi_in_parallel.sh
```

In this step we do all the inference with dadi. Detailed comments can be found within in the Python scripts in this bash script.  We've also included the contents of the domesticationDFE/results/data_frames/ folder after running the dadi and polyDFE inferences.  We've made sure to include all the data frames you'll need to reproduce all the figures in the manuscript (except for Supplementary Figures 6 and 7, which require the contents of the domesticationDFE/results/VCF/large folder).

# Stopped HERE




## Usage 





