import dadi
import sys
import matplotlib.pyplot as plt
import csv
import pickle
import dadi.DFE as DFE
import domestication_new_dadi_functions as new_models
demography_domestication = new_models.Domestication_demography
demography_complex_domestication = new_models.Domestication_flexible_demography


migration = sys.argv[1]
possel = sys.argv[2]
change = sys.argv[3]
# iteration  = sys.argv[4]
site = "syn"

# You can visualize the comparison between the data and model

# Using demographic inference information to generate the model 
with open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'.2d.flexible_dem_'+site+'_best_fit.txt', newline = '') as best_demo:  
    best_demo_reader = csv.reader(best_demo, delimiter='\t')
    for popt in best_demo_reader:
        print(popt)

theta_syn = float(popt[13]) 
popt = popt[1:13]
popt = list(map(float, popt))
print(theta_syn)
print('Best-fit parameters: {0}'.format(popt))

# Loading data and setting up what we need for the model  #sim_0-200-0.25-nsyn.2d.fs
data_fs = dadi.Spectrum.from_file('../results/SFS/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_subsampled.2d.fs')
ns = data_fs.sample_sizes
print(ns)
pts_l = [max(ns)+20, max(ns)+30, max(ns)+40]

func = demography_complex_domestication
func = dadi.Numerics.make_anc_state_misid_func(func)
func_ex = dadi.Numerics.make_extrap_log_func(func)

# Syn SFS plot
model_fs = func_ex(popt, ns, pts_l)

import pylab
pylab.figure(figsize=(8,6))
dadi.Plotting.plot_2d_comp_multinom(model_fs, data_fs, pop_ids =('Wild','Domesticated'), show=False)

pylab.savefig('../plots/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_data_vs_flex_model.2d.png', dpi=250)

