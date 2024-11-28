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

