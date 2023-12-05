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
# iteration  = sys.argv[4]
site = "syn"

# Load the data to just get ns
data = dadi.Spectrum.from_file('../results/SFS/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_subsampled.2d.fs')
ns = data.sample_sizes
print(ns)
pts_l = [max(ns)+440, max(ns)+450, max(ns)+460]
add_gammas, bounds, gamma_pts = [0,1,10,100], (1e-4, 2000), 100

# with open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'.2d.complex_dem_'+site+'_best_fit.txt', newline = '') as best_demo: 
with open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'.2d.flexible_dem_'+site+'_best_fit.txt', newline = '') as best_demo:   
    best_demo_reader = csv.reader(best_demo, delimiter='\t')
    for popt in best_demo_reader:
        print(popt)

# theta_syn = float(popt[18]) 
# popt = popt[1:16]
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
# fid = open('s1_sfs{0}.bpkl'.format(filenames[scenario_ii]), 'wb')
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
# fid = open('s2_sfs{0}.bpkl'.format(filenames[scenario_ii]), 'wb')
pickle.dump(s2, fid, protocol=2)
fid.close()






# def Domestication(params, ns, pts):
    
#     nu1,nuPre,nu2,TPre,T,T1,m,gamma1,gamma2 = params
   
#     xx = Numerics.default_grid(pts)
   
#     phi = PhiManip.phi_1D(xx, gamma=gamma1)
#     phi = Integration.one_pop(phi, xx, TPre, nu=nu1, gamma=gamma1)
#     phi = PhiManip.phi_1D_to_2D(xx, phi)
   
#     phi = Integration.two_pops(phi, xx, T, nu1, nuPre, m21=m, gamma1=gamma1, gamma2=gamma2)
#     phi = Integration.two_pops(phi, xx, T1, nu1, nu2,  m21=m, gamma1=gamma1, gamma2=gamma2)
#     fs = Spectrum.from_phi(phi, ns, (xx,xx))
#     return fs
   
   
# def Domestication_one(params, ns, pts):
#     newpars = list(params[:-1]) + [params[-1], params[-1]]
#     return Domestication(newpars, ns, pts)

# if __name__ == '__main__':
# 	popts = {1:[0.9882964648680578,	0.07599335011848302,	1.3973755767634704,	94.5890220722061,	8.616984305218493,	0.0008301767751998229,	0.13937775449604375],
# 	2:[0.499615716	,0.280029145,	3.240488823,	15	,4.713920005	,0.115240363	,0.006195742],
# 	3:[0.6997682218672089,0.10617519634660136,1.8765497536517077, 6.726575918665505,17.9303021134617,	0.23720143822954767	,0.039113528770318436],
# 	4:[1.2746126700706815,	0.06884699304154243	,1.3225397657188886,	0.0010099848371629288	,6.238253507517234,	0.014968766126308514	,0.2893266527238452],
# 	5:[0.8558178450700078,	0.07503547574906316,	0.886956859470374,	0.001010628107083513,	9.807128539408806	,0.03106850646103766	,0.1848006089791323],
# 	6:[0.62229691533924,	0.10024467257674997	,0.6898718139331717,	0.0010194475384830402	,6.260742923328892,	0.04734610390332232,	0.13167270824123023],
# 	7:[0.651392469,	0.0100809, 2.307906263,	0.001002079,	13.46668305,	0.174889179,	0.005734094],
# 	8:[0.8279576047003072,	0.010000000000000004,	3.3676267996388596,	0.0010000000000000002	,7.959480077999102,	0.2106674564333588,	0.006877128400215741],
# 	9:[0.584616797,	0.010077501,	1.858109455,	0.001011808	, 16.51722028	,0.143000171	,0.004771849],
# 	10:[0.640732339,	0.010360706,	3.555215676,	0.001004896,	7.603493467	, 0.174441889,	0.006071097]}
	
# 	# Rearrange popts as Ioanna did
# 	for scen, popt in popts.items():
# 	    popts[scen]=[popt[0],popt[1],popt[2],popt[4],popt[5],popt[6],popt[3]]
	
# 	ns = [40, 40]
# 	pts_l = [max(ns)+140, max(ns)+150, max(ns)+160]
# 	add_gammas, bounds, gamma_pts = [0,1,10], (1e-4, 2000), 50

# 	## For quick testing
# 	# bounds, gamma_pts, add_gammas, pts_l = (1e-4, 2), 2, [1], [41,42,43]

# 	# Preserve Ioanna's filename convention
# 	filenames = {1:19, 10:86, 2:28, 3:98, 4:11, 5:31, 6:64, 7:16, 8:53, 9:33}
	   
# 	s1 = Cache1D(popts[scenario_ii], ns, Domestication_one, pts=pts_l, gamma_pts=gamma_pts,
# 	             gamma_bounds=bounds, additional_gammas=add_gammas, verbose=True, mp=True)

# 	fid = open('s1_sfs{0}.bpkl'.format(filenames[scenario_ii]), 'wb')
# 	pickle.dump(s1, fid, protocol=2)
# 	fid.close()
	   
# 	s2 = Cache2D(popts[scenario_ii], ns, Domestication, pts=pts_l, gamma_pts=gamma_pts,
# 	                  gamma_bounds=bounds, additional_gammas=add_gammas, verbose=True, mp=True)
	  
# 	fid = open('s2_sfs{0}.bpkl'.format(filenames[scenario_ii]), 'wb')
# 	pickle.dump(s2, fid, protocol=2)
# 	fid.close()
	
