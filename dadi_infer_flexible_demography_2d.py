import dadi
import nlopt
import sys
import random
import domestication_new_dadi_functions as new_models
demography_domestication = new_models.Domestication_flexible_demography

random.seed(12345)

migration = sys.argv[1]
possel = sys.argv[2]
change = sys.argv[3]
# iteration  = sys.argv[4]
site = "syn"


# Load the data
data = dadi.Spectrum.from_file('../results/SFS/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_subsampled.2d.fs')
ns = data.sample_sizes

# These are the grid point settings will use for extrapolation.
pts_l = [50,60,70]

# The Demographics1D and Demographics2D modules contain a few simple models,
# mostly as examples. We could use one of those.
# func = dadi.Demographics2D.split_mig
# Instead, we'll work with our custom model
func = demography_domestication

# Now let's optimize parameters for this model.

# The upper_bound and lower_bound lists are for use in optimization.
# Occasionally the optimizer will try wacky parameter values. We in particular
# want to exclude values with very long times, very small population sizes, or
# very high migration rates, as they will take a long time to evaluate.
# Parameters are: Tpre, nuPre, Tdiv, nu1div, nu2div, T1F, T2F, nu1F, nu2F, mw2d, md2w = params

upper_bounds =    [1,   5,     1,    2,      1,      1,   2,   10,   10,   1,    1,   0.5]
lower_bounds =    [0,   1e-4,  0,    1e-4,   1e-5,   0,   0,   1e-4, 1e-4, 0,    0,   0]

# This is our initial guess for the parameters, which is based on the average best values for scenarios with very weak and weak positive selection.
params =   [0.1812, 1.6429, 0.0522, 0.6668, 0.0730, 0.2267, 0.2845, 2.1584, 1.8396, 0.0283, 0.0411, 0.0035]
fix_params = [None, None, None, None, None, None, None, None, None, None, None, None]

# Make the extrapolating version of our demographic model function.
func = dadi.Numerics.make_anc_state_misid_func(func)
func_ex = dadi.Numerics.make_extrap_log_func(func)


try:
    fid = open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'.2d.flexible_dem_syn_fits.txt','a')
except:
    fid = open('../results/dadi_outputs/sim_'+migration+'-'+possel+'-'+change+'.2d.flexible_dem_syn_fits.txt','w')


for i in range(100):
    p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bounds,lower_bound=lower_bounds)

    # Run optimization
    # At the end of the optimization you will get the
    # optimal parameters and log-likelihood.
    # You can modify verbose to watch how the optimizer behaves,
    # what number you pass it how many evaluations are done
    # before the evaluation is printed.
    popt, ll_model = dadi.Inference.opt(p0, data, func_ex, pts_l, lower_bound=lower_bounds, upper_bound=upper_bounds, algorithm=nlopt.LN_BOBYQA, maxeval=600, fixed_params=fix_params, verbose=True)

    # Find the synonymous theta
    model_fs = func_ex(popt, ns, pts_l)
    theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, data)

    # Optional save method 1:
    # Write results to fid
    res = [ll_model] + list(popt) + [theta0]
    fid.write('\t'.join([str(ele) for ele in res])+'\n')


fid.close()


# # Perturb our parameters before optimization. This does so by taking each
# # parameter a up to a factor of two up or down.
# p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,
#                               lower_bound=lower_bound)
# # Do the optimization. By default we assume that theta is a free parameter,
# # since it's trivial to find given the other parameters. If you want to fix
# # theta, add a multinom=False to the call.
# # The maxiter argument restricts how long the optimizer will run. For real 
# # runs, you will want to set this value higher (at least 10), to encourage
# # better convergence. You will also want to run optimization several times
# # using multiple sets of intial parameters, to be confident you've actually
# # found the true maximum likelihood parameters.
# print('Beginning optimization ************************************************')
# popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l, 
#                                    lower_bound=lower_bound,
#                                    upper_bound=upper_bound,
#                                    verbose=len(p0), maxiter=3)
# # The verbose argument controls how often progress of the optimizer should be
# # printed. It's useful to keep track of optimization process.
# print('Finshed optimization **************************************************')

# # These are the actual best-fit model parameters, which we found through
# # longer optimizations and confirmed by running multiple optimizations.
# # We'll work with them through the rest of this script.


# popt = [1.880, 0.0724, 1.764, 0.930, 0.363, 0.112]
# print('Best-fit parameters: {0}'.format(popt))

# # Calculate the best-fit model AFS.
# model = func_ex(popt, ns, pts_l)
# # Likelihood of the data given the model AFS.
# ll_model = dadi.Inference.ll_multinom(model, data)
# print('Maximum log composite likelihood: {0}'.format(ll_model))
# # The optimal value of theta given the model.
# theta = dadi.Inference.optimal_sfs_scaling(model, data)
# print('Optimal value of theta: {0}'.format(theta))

# # Plot a comparison of the resulting fs with the data.
# import pylab
# pylab.figure(figsize=(8,6))
# dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3, pop_ids =('Domesticated','Wild'), show=False)
# # Save the figure
# pylab.savefig('../plots/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+iteration+'-'+site+'_data_vs_model.2d.png', dpi=250)
