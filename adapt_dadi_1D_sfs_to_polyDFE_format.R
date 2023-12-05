library("data.table")
library("dplyr")

args = commandArgs(trailingOnly=TRUE)

#args[1] = "0" 
#args[2] = "2"
#args[3] = "0"
#args[4] = "1"

MIGRATION = args[1]
POSSEL = args[2]
CHANGE = args[3]
ITER = args[4]
header = "1 1 20"

for (POP in c("wild", "domesticated"))
{
	# POP = "wild"
	SFS_syn   = fread(paste0("../results/SFS/by_dadi/Bootstraps/sim_", MIGRATION, "-", POSSEL, "-", CHANGE, "-syn_boot.",  ITER, ".", POP, "_subsampled.1d.fs" ), skip=1, nrows=1, header = F)
	SFS_nsyn  = fread(paste0("../results/SFS/by_dadi/Bootstraps/sim_", MIGRATION, "-", POSSEL, "-", CHANGE, "-nsyn_boot.", ITER, ".", POP, "_subsampled.1d.fs" ), skip=1, nrows=1, header = F)

	SFS_syn = dplyr::select(SFS_syn, -V1)  
	SFS_nsyn = dplyr::select(SFS_nsyn, -V1)  

	SFS_syn = append(SFS_syn,  400000*100, after = 19) #400,000 is the total number of synonymous "sites" in our slim simulations, since we perform 100 replicates we multiply both numbers. This must be changed if required
	SFS_syn = append(SFS_syn,  400000*100, after = 21)

	SFS_nsyn = append(SFS_nsyn,  800000*100, after = 19) #800,000 is the total number of nonsynonymous "sites" in our slim simulations, since we perform 100 replicates we multiply both numbers. This must be changed if required
	SFS_nsyn = append(SFS_nsyn,  800000*100, after = 21)

	SFS_syn = unlist(SFS_syn)
	SFS_nsyn = unlist(SFS_nsyn)

	print(SFS_syn)
	print(SFS_nsyn)

	write(header,   paste0("../results/SFS/for_polyDFE/sim_", MIGRATION, "-", POSSEL, "-", CHANGE, "-", ITER, "-", POP, ".fs"), sep="\t", append=F)
	write(SFS_syn,  paste0("../results/SFS/for_polyDFE/sim_", MIGRATION, "-", POSSEL, "-", CHANGE, "-", ITER, "-", POP, ".fs"), sep="\t", append=T, ncolumns = 22)
	write(SFS_nsyn, paste0("../results/SFS/for_polyDFE/sim_", MIGRATION, "-", POSSEL, "-", CHANGE, "-", ITER, "-", POP, ".fs"), sep="\t", append=T, ncolumns = 22)
}