import pickle
import nlopt
import sys
import random
import csv
from dadi.DFE import *
from dadi import Numerics, Integration, PhiManip, Spectrum
from domestication_new_dadi_functions import *

random.seed(12345)

migration = sys.argv[1]
possel = sys.argv[2]
change = sys.argv[3]
site = "syn"

# Load the data to just get ns
data = dadi.Spectrum.from_file('../results/SFS/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_subsampled.2d.fs')
ns = data.sample_sizes
print(ns)
pts_l = [max(ns)+440, max(ns)+450, max(ns)+460]
add_gammas, bounds, gamma_pts = [0,1,10,100], (1e-4, 2000), 100

with open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'.2d.flexible_dem_'+site+'_best_fit.txt', newline = '') as best_demo:   
    best_demo_reader = csv.reader(best_demo, delimiter='\t')
    for popt in best_demo_reader:
        print(popt)

theta_syn = float(popt[13]) 
popt = popt[1:12]
popt = list(map(float, popt))
print('Best-fit parameters: {0}'.format(popt))


s1 = Cache1D(popt, ns, Domestication_one_pop, pts=pts_l, gamma_pts=gamma_pts, gamma_bounds=bounds, additional_gammas=add_gammas, verbose=True, cpus=None)

# Check if the cached spectra have any large negative values
if (s1.spectra<0).sum() > 0:
    print(
        '!!!WARNING!!!\nPotentially large negative values!\nMost negative value is: '+str(s1.spectra.min())+
        '\nIf negative values are very negative (<-0.001), rerun with larger values for pts_l'
        )

# Save the cache with pickle
fid = open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'_1d_cache.bpkl', 'wb')
pickle.dump(s1, fid, protocol=2)
fid.close()
   
s2 = Cache2D(popt, ns, Domestication_flexible_demography_cache, pts=pts_l, gamma_pts=gamma_pts, gamma_bounds=bounds, additional_gammas=add_gammas, verbose=True, cpus=None)

# Check if the cached spectra have any large negative values
if (s2.spectra<0).sum() > 0:
    print(
        '!!!WARNING!!!\nPotentially large negative values!\nMost negative value is: '+str(s2.spectra.min())+
        '\nIf negative values are very negative (<-0.001), rerun with larger values for pts_l'
        )

# Save the cache with pickle  
fid = open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'_2d_cache.bpkl', 'wb')
pickle.dump(s2, fid, protocol=2)
fid.close()
