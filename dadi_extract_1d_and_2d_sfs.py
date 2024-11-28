import dadi
import nlopt
import sys
import random

random.seed(12345)

migration = sys.argv[1]
possel = sys.argv[2]
change = sys.argv[3]
site = sys.argv[4]


### 2D-SFS ##
pop_ids, ns = ['Wild','Domesticated'], [40,40]

datafile = '../results/VCF/sim_'+migration+'-'+possel+'-'+change+'_'+site+'.one_head.vcf'

dd = dadi.Misc.make_data_dict_vcf(datafile, 'Wild_Domesticated_individuals.txt')

fs = dadi.Spectrum.from_data_dict(dd, pop_ids, ns, polarized=True)

fs.to_file('../results/SFS/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+site+'.2d.fs')

fs_folded = dadi.Spectrum.from_data_dict(dd, pop_ids, ns, polarized=False)

new_ns = [20, 20]
dd_subsample = dadi.Misc.make_data_dict_vcf(datafile, 'Wild_Domesticated_individuals.txt', subsample={'Domesticated': ns[0]//2, 'Wild': ns[1]//2})
fs_subsample = dadi.Spectrum.from_data_dict(dd_subsample, pop_ids, new_ns)
fs_subsample.to_file('../results/SFS/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_subsampled.2d.fs')


# Bootstrap parameters #
Nboot, chunk_size = 100, 1.2e6 #this is the size of the simulated genomic region. So it will take each independent slim replicate as a chunk.
chunks = dadi.Misc.fragment_data_dict(dd_subsample, chunk_size)
boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, Nboot, pop_ids, new_ns)

# Print boostrapped FS #
for i in range(len(boots)):
    boots[i].to_file('../results/SFS/by_dadi/Bootstraps/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_boot.{0}.subsampled.2d.fs'.format(str(i)))



import matplotlib.pyplot as plt

fig = plt.figure(1, figsize=(12,10))
fig.clear()

# Note that projection creates fractional entries in the spectrum,
# so vmin < 1 is sensible.
ax = fig.add_subplot(2,2,1)
dadi.Plotting.plot_single_2d_sfs(fs, vmin=1e-2, ax=ax)
ax.set_title('Orignal data')

ax = fig.add_subplot(2,2,2)
dadi.Plotting.plot_single_2d_sfs(fs_folded, vmin=1e-2, ax=ax)
ax.set_title('Folded original data')

ax = fig.add_subplot(2,2,3)
dadi.Plotting.plot_single_2d_sfs(fs_subsample, vmin=1e-2, ax=ax)
ax.set_title('Subsampled original data')

ax = fig.add_subplot(2,2,4)
dadi.Plotting.plot_single_2d_sfs(boots[0], vmin=1e-2, ax=ax)
ax.set_title('Bootstrap subsampled data')

fig.tight_layout()
# plt.show()
fig.savefig('../plots/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+site+'.2d.png')



### 1D-SFS - Wild ##
pop_ids, ns = ['Wild'], [40]

datafile = '../results/VCF/sim_'+migration+'-'+possel+'-'+change+'_'+site+'.one_head.vcf'

dd = dadi.Misc.make_data_dict_vcf(datafile, 'Wild_Domesticated_individuals.txt')

fs = dadi.Spectrum.from_data_dict(dd, pop_ids, ns, polarized=True)

fs.to_file('../results/SFS/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_wild.1d.fs')

fs_folded = dadi.Spectrum.from_data_dict(dd, pop_ids, ns, polarized=False)

new_ns = [20]
dd_subsample = dadi.Misc.make_data_dict_vcf(datafile, 'Wild_Domesticated_individuals.txt')
fs_subsample = dadi.Spectrum.from_data_dict(dd_subsample, pop_ids, new_ns)
fs_subsample.to_file('../results/SFS/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_wild_subsampled.1d.fs')


# Bootstrap parameters #
Nboot, chunk_size = 100, 1.2e6 #this is the size of the simulated genomic region. So it will take each independent slim replicate as a chunk.
chunks = dadi.Misc.fragment_data_dict(dd_subsample, chunk_size)
boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, Nboot, pop_ids, new_ns)

# Print boostrapped FS #
for i in range(len(boots)):
    boots[i].to_file('../results/SFS/by_dadi/Bootstraps/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_boot.{0}.wild_subsampled.1d.fs'.format(str(i)))



import matplotlib.pyplot as plt

fig = plt.figure(1, figsize=(10,6))
fig.clear()

# Note that projection creates fractional entries in the spectrum,
# so vmin < 1 is sensible.
#ax = fig.add_subplot(2,2,1)
#dadi.Plotting.plot_1d_fs(fs, fig_num=1)
#ax.set_title('Orignal data')

#ax = fig.add_subplot(2,2,2)
#dadi.Plotting.plot_1d_fs(fs_folded, fig_num=2)
#ax.set_title('Folded original data')

#ax = fig.add_subplot(2,2,3)
dadi.Plotting.plot_1d_fs(fs_subsample)
ax.set_title('Subsampled original data')

#ax = fig.add_subplot(2,2,4)
#dadi.Plotting.plot_1d_fs(boots[0], fig_num=4)
#ax.set_title('Bootstrap subsampled data')

fig.tight_layout()
# plt.show()
fig.savefig('../plots/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+site+'.wild_1d.png')



### 1D-SFS - Domesticated ##
pop_ids, ns = ['Domesticated'], [40]

datafile = '../results/VCF/sim_'+migration+'-'+possel+'-'+change+'_'+site+'.one_head.vcf'

dd = dadi.Misc.make_data_dict_vcf(datafile, 'Wild_Domesticated_individuals.txt')

fs = dadi.Spectrum.from_data_dict(dd, pop_ids, ns, polarized=True)

fs.to_file('../results/SFS/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_domesticated.1d.fs')

fs_folded = dadi.Spectrum.from_data_dict(dd, pop_ids, ns, polarized=False)

new_ns = [20]
dd_subsample = dadi.Misc.make_data_dict_vcf(datafile, 'Wild_Domesticated_individuals.txt')
fs_subsample = dadi.Spectrum.from_data_dict(dd_subsample, pop_ids, new_ns)
fs_subsample.to_file('../results/SFS/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_domesticated_subsampled.1d.fs')


# Bootstrap parameters #
Nboot, chunk_size = 100, 1.2e6 #this is the size of the simulated genomic region. So it will take each independent slim replicate as a chunk.
chunks = dadi.Misc.fragment_data_dict(dd_subsample, chunk_size)
boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, Nboot, pop_ids, new_ns)

# Print boostrapped FS #
for i in range(len(boots)):
    boots[i].to_file('../results/SFS/by_dadi/Bootstraps/sim_'+migration+'-'+possel+'-'+change+'-'+site+'_boot.{0}.domesticated_subsampled.1d.fs'.format(str(i)))



import matplotlib.pyplot as plt

fig = plt.figure(1, figsize=(10,6))
fig.clear()

# Note that projection creates fractional entries in the spectrum,
# so vmin < 1 is sensible.
#ax = fig.add_subplot(2,2,1)
#dadi.Plotting.plot_1d_fs(fs, fig_num=1)
#ax.set_title('Orignal data')

#ax = fig.add_subplot(2,2,2)
#dadi.Plotting.plot_1d_fs(fs_folded, fig_num=2)
#ax.set_title('Folded original data')

#ax = fig.add_subplot(2,2,3)
dadi.Plotting.plot_1d_fs(fs_subsample)
ax.set_title('Subsampled original data')

#ax = fig.add_subplot(2,2,4)
#dadi.Plotting.plot_1d_fs(boots[0], fig_num=4)
#ax.set_title('Bootstrap subsampled data')

fig.tight_layout()
# plt.show()
fig.savefig('../plots/by_dadi/sim_'+migration+'-'+possel+'-'+change+'-'+site+'.domesticated_1d.png')


