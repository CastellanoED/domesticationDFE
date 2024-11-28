### dadi_uncertainty 
# Generate 100 bootstrap datasets, by dividing the genome into 2 Mb chunks and
# resampling across those chunks.
import dadi
import random
import nlopt
import sys
import csv
import domestication_new_dadi_functions as new_models
demography_domestication = new_models.Domestication_flexible_demography


random.seed(12345)

migration = sys.argv[1]
possel = sys.argv[2]
change = sys.argv[3]
site = "syn"

pop_ids, ns = ['Wild','Domesticated'], [40,40]

datafile = '../results/VCF/sim_'+migration+'-'+possel+'-'+change+'_'+site+'.one_head.vcf'
dd = dadi.Misc.make_data_dict_vcf(datafile, 'Wild_Domesticated_individuals.txt')

new_ns = [20, 20]
dd_subsample = dadi.Misc.make_data_dict_vcf(datafile, 'Wild_Domesticated_individuals.txt', subsample={'Domesticated': ns[0]//2, 'Wild': ns[1]//2})
fs = dadi.Spectrum.from_data_dict(dd_subsample, pop_ids, new_ns)


# Bootstrap parameters #
Nboot, chunk_size = 100, 1.2e6 #this is the size of the simulated genomic region. So it will take each independent slim replicate as a chunk.
chunks = dadi.Misc.fragment_data_dict(dd, chunk_size)
boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, Nboot, pop_ids, new_ns)


# Define the demographic model
# Single population demographic models are in dadi.Demographics1D
# Two population demographic models are in dadi.Demographics2D
# demo_model = dadi.Demographics1D.two_epoch
demo_model = demography_domestication
# demo_model = dadi.Demographics1D.three_epoch
# demo_model = three_epoch_inbr


# If the data is unfolded (the ancestral allele was known), as the example data is
# Wrap the demographic model in a function that adds a parameter to estimate the 
# rate of misidentification.
demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)


# Wrap the demographic model in a function that utilizes grid points
# which increases dadi's ability to more accurately generate a model
# frequency spectrum.
demo_model_ex = dadi.Numerics.make_extrap_log_func(demo_model)


# Define the grid points based on the sample size
# For smaller data (largest sample size ~ <100) [ns+20, ns+30, ns+40] is a good starting point
# For larger data (largest sample size ~ >=100) [ns+100, ns+110, ns+120] is a good starting point
# pts_l = [max(ns)+120, max(ns)+130, max(ns)+140]
pts_l = [max(ns)+20, max(ns)+30, max(ns)+40]


# Define the bestfit parameters
with open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'.2d.flexible_dem_'+site+'_best_fit.txt', newline = '') as best_demo:   
    best_demo_reader = csv.reader(best_demo, delimiter='\t')
    for popt_demo in best_demo_reader:
        print(popt_demo)

theta_syn = float(popt_demo[13]) 
popt_demo = popt_demo[1:13]
popt_demo = list(map(float, popt_demo))
print('Best-fit parameters: {0}'.format(popt_demo))


# Godambe uncertainties
# Will contain uncertainties for the
# estimated demographic parameters and theta.

# Start a file to contain the confidence intervals
fi = open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'_demographic_confidence_intervals.txt','wt')
# fi.write('Optimized parameters:{0}\n\n'.format(popt_demo))

# we want to try a few different step sizes (eps) to see if
# uncertainties vary wildly with changes to step size.
for eps in [0.01, 0.001, 0.0001]:
    uncerts_adj_demo= dadi.Godambe.GIM_uncert(demo_model_ex, pts_l, boots, popt_demo, fs, eps=eps)
    # fi.write('Estimated 95% uncerts (with step size '+str(eps)+'):{0}\n'.format(1.96*uncerts_adj_demo[:-1]))
    fi.write('Lower_bound '+str(eps)+' {0}\n'.format(popt_demo-1.96*uncerts_adj_demo[:-1]))
    fi.write('Upper_bound '+str(eps)+' {0}\n'.format(popt_demo+1.96*uncerts_adj_demo[:-1]))
fi.close()



