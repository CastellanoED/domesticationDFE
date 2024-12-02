# Dadi workflow example 2 - Infer the demographic model
import dadi
import sys
import matplotlib.pyplot as plt
import csv
import pickle
import dadi.DFE as DFE
import domestication_new_dadi_functions as new_models
demography_domestication = new_models.Domestication_demography
demography_complex_domestication = new_models.Vourlaki_mixture

# Make a variable to store the name of the dataset you are working with
# so that you can esaily change it to work on different datasets
# dataset = 'ACG_T'

migration = sys.argv[1]
possel = sys.argv[2]
change = sys.argv[3]
sb = sys.argv[4]

# iteration  = sys.argv[4]
site = "nsyn"

# You can visualize the comparison between the data and model

# Using demographic inference information to generate the model

# sim_0.01-200-0.25-nsyn_fits_with_Sb=0_two_pb.txt
with open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_with_Sb='+sb+'_best_fit.txt', newline = '') as best_demo:  
    best_demo_reader = csv.reader(best_demo, delimiter='\t')
    for popt in best_demo_reader:
        print(popt)

theta_nsyn = float(popt[8]) 
popt = popt[1:8]
popt = list(map(float, popt))
print(theta_nsyn)
print('Best-fit parameters: {0}'.format(popt))

# Loading data and setting up what we need for the model  #sim_0-200-0.25-nsyn.2d.fs
data_fs = dadi.Spectrum.from_file('../results/SFS/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_subsampled.2d.fs')
ns = data_fs.sample_sizes
print(ns)
pts_l = [max(ns)+20, max(ns)+30, max(ns)+40]

s1 = pickle.load(open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'_1d_cache.bpkl','rb'))
s2 = pickle.load(open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'_2d_cache.bpkl','rb'))

func = demography_complex_domestication
func = dadi.Numerics.make_anc_state_misid_func(func)

# NSyn SFS plot
model_fs = func(popt, ns, s1, s2, theta_nsyn, pts_l)

import pylab
pylab.figure(figsize=(8,6))
dadi.Plotting.plot_2d_comp_multinom(model_fs, data_fs, pop_ids =('Wild','Domesticated'), show=False)

pylab.savefig('../plots/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+sb+'-'+site+'_data_vs_vourlaki_model.2d.png', dpi=250)

