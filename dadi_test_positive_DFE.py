# Dadi workflow example 4 - Infer the DFE
import dadi
import dadi.DFE as DFE
import pickle
import nlopt
import sys
import csv
import domestication_new_dadi_functions as new_models
dfe_domestication = new_models.Vourlaki_mixture
dfe_domestication_uniq_ppos = new_models.Vourlaki_mixture_uniq_ppos

migration = sys.argv[1]
possel = sys.argv[2]
change = sys.argv[3]
# iteration  = sys.argv[4]
site = "nsyn"

# Load the original data                                                     
data = dadi.Spectrum.from_file('../results/SFS/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_subsampled.2d.fs')
ns = data.sample_sizes

### TEST 4 DIFFERENT STRENGTHS OF POSITIVE SELECTION and two pb
for new_possel in [2,20,200,0]: 
    print(possel)
    print(new_possel)
    # Using inference information to generate the model
    with open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'.2d.flexible_dem_syn_best_fit.txt', newline = '') as best_demo:   
        best_demo_reader = csv.reader(best_demo, delimiter='\t')
        for popt in best_demo_reader:
            print(popt)

    theta_syn = float(popt[13]) 
    nsyn_syn_ratio = float(popt[14])
    popt = popt[1:12]
    popt = list(map(float, popt))
    print('Best-fit parameters: {0}'.format(popt))

    gamma_pos = int(new_possel)/2
    if gamma_pos > 0: 
        ppos_wild = 0.1/gamma_pos
    if gamma_pos == 0: 
        ppos_wild = 0

    # Calculate the theta for the nonsynonymous data based
    # on the ratio of nonsynonymous mutations to synonymous mutations.
    # You will probably be get this ratio from the paper or calculate it,
    # but typically it is larger grater than 1.
    theta = theta_syn * nsyn_syn_ratio

    print(gamma_pos)
    print(ppos_wild)

    print(theta_syn)
    print(nsyn_syn_ratio)
    print(theta)

    # Load the cache of spectra      #sim_0-200-0_1d_cache.bpkl
    s1 = pickle.load(open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'_1d_cache.bpkl','rb'))
    s2 = pickle.load(open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'_2d_cache.bpkl','rb'))

    # Optimization for the DFE requires extra arguments that
    # are not included in the optimizer function, so we need
    # to define them ourselves.
    func_args = [s1, s2, theta]

    # Choose starting parameters for inference
    params = [0.3, 100, ppos_wild, gamma_pos, float(change), ppos_wild, 0.01]

    lower_bound=[0.01,0.01, 0,  0,0,0,0]
    upper_bound=[10,  2000, 1,100,1,1,1]

    if gamma_pos > 0: #only the strength is fix, the probability can vary
        fixed_params=[None,None,None,gamma_pos,None,None,None]
    #params =  alpha, beta, ppos_wild, gamma_pos, pchange, pchange_pos 

    if gamma_pos == 0: #no postive DFE
        fixed_params=[None,None,0,0,None,0,None]


    func_ex = dadi.Numerics.make_anc_state_misid_func(dfe_domestication)


    # Optional: create a file to store fits in
    # For running on the HPC, it is a good idea to
    # check if file is made before making it so
    # that you don't overwrite other results
    try:
    	fid = open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_fits_with_Sb='+str(new_possel)+'_two_pb.txt','a')
    except:
    	fid = open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_fits_with_Sb='+str(new_possel)+'_two_pb.txt','w')


    # Perturb parameters
    # Optimizers dadi uses are mostly deterministic
    # so you will want to randomize parameters for each optimization.
    # It is recommended to optimize at least 20 time if you only have
    # your local machin, 100 times if you have access to an HPC.
    # If you want a single script to do multiple runs, you will want to
    # start a for loop here
    for i in range(100):
        p0 = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

        # Run optimization
        # At the end of the optimization you will get the
        # optimal parameters and log-likelihood.
        # You can modify verbose to watch how the optimizer behaves,
        # what number you pass it how many evaluations are done
        # before the evaluation is printed.
        # For the DFE, because we calculated the theta_nsyn, we want to
        # set multinom=False.

        popt, ll = dadi.Inference.opt(p0, data, func_ex, pts=None, func_args=func_args, lower_bound=lower_bound, upper_bound=upper_bound, fixed_params=fixed_params, multinom=False, maxeval=1000, verbose=True)

        # Optional save method 1:
        # Write results to fid
        res = [ll] + list(popt) + [theta]

        fid.write('\t'.join([str(ele) for ele in res])+'\n')
        # fid.write(os.lines ep)

    # Optional save method 1: 
    # Close the file
    fid.close()


### TEST 4 DIFFERENT STRENGTHS OF POSITIVE SELECTION and one pb
for new_possel in [2,20,200,0]: 
    print(possel)
    print(new_possel)
    # Using inference information to generate the model
    with open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'.2d.flexible_dem_syn_best_fit.txt', newline = '') as best_demo:   
        best_demo_reader = csv.reader(best_demo, delimiter='\t')
        for popt in best_demo_reader:
            print(popt)

    theta_syn = float(popt[13]) 
    nsyn_syn_ratio = float(popt[14])
    popt = popt[1:12]
    popt = list(map(float, popt))
    print('Best-fit parameters: {0}'.format(popt))

    gamma_pos = int(new_possel)/2
    if gamma_pos > 0: 
        ppos_wild = 0.1/gamma_pos
    if gamma_pos == 0: 
        ppos_wild = 0

    # Calculate the theta for the nonsynonymous data based
    # on the ratio of nonsynonymous mutations to synonymous mutations.
    # You will probably be get this ratio from the paper or calculate it,
    # but typically it is larger grater than 1.
    theta = theta_syn * nsyn_syn_ratio

    print(gamma_pos)
    print(ppos_wild)

    print(theta_syn)
    print(nsyn_syn_ratio)
    print(theta)

    # Load the cache of spectra      #sim_0-200-0_1d_cache.bpkl
    s1 = pickle.load(open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'_1d_cache.bpkl','rb'))
    s2 = pickle.load(open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'_2d_cache.bpkl','rb'))

    # Optimization for the DFE requires extra arguments that
    # are not included in the optimizer function, so we need
    # to define them ourselves.
    func_args = [s1, s2, theta]

    # Choose starting parameters for inference
    params = [0.3, 100, ppos_wild, gamma_pos, float(change), 0.01]

    lower_bound=[0.01,0.01, 0,  0,0,0]
    upper_bound=[10,  2000, 1,100,1,1]

    if gamma_pos > 0: #only the strength is fix, the probability can vary
        fixed_params=[None,None,None,gamma_pos,None,None]
    #params =  alpha, beta, ppos_wild, gamma_pos, pchange 

    if gamma_pos == 0: #no postive DFE
        fixed_params=[None,None,0,0,None,None]

    func_ex = dadi.Numerics.make_anc_state_misid_func(dfe_domestication_uniq_ppos)


    # Optional: create a file to store fits in
    # For running on the HPC, it is a good idea to
    # check if file is made before making it so
    # that you don't overwrite other results
    try:
        fid = open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_fits_with_Sb='+str(new_possel)+'_uniq_pb.txt','a')
    except:
        fid = open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_fits_with_Sb='+str(new_possel)+'_uniq_pb.txt','w')


    # Perturb parameters
    # Optimizers dadi uses are mostly deterministic
    # so you will want to randomize parameters for each optimization.
    # It is recommended to optimize at least 20 time if you only have
    # your local machin, 100 times if you have access to an HPC.
    # If you want a single script to do multiple runs, you will want to
    # start a for loop here
    for i in range(100):
        p0 = dadi.Misc.perturb_params(params, fold=2, upper_bound=upper_bound, lower_bound=lower_bound)

        # Run optimization
        # At the end of the optimization you will get the
        # optimal parameters and log-likelihood.
        # You can modify verbose to watch how the optimizer behaves,
        # what number you pass it how many evaluations are done
        # before the evaluation is printed.
        # For the DFE, because we calculated the theta_nsyn, we want to
        # set multinom=False.

        popt, ll = dadi.Inference.opt(p0, data, func_ex, pts=None, func_args=func_args, lower_bound=lower_bound, upper_bound=upper_bound, fixed_params=fixed_params, multinom=False, maxeval=1000, verbose=True)

        # Optional save method 1:
        # Write results to fid
        res = [ll] + list(popt) + [theta]

        fid.write('\t'.join([str(ele) for ele in res])+'\n')
        # fid.write(os.lines ep)

    # Optional save method 1: 
    # Close the file
    fid.close()