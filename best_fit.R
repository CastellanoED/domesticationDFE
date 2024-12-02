library("data.table")
library("stringr")
args = commandArgs(trailingOnly=TRUE)

# args[1] = 0
# args[2] = 2
# args[3] = 0
# args[4] = "true" #or complex


if (args[4] == "true")    
{ 
	data = read.delim(paste0("../results/dadi_outputs/sim_", args[1], "-", args[2], "-", args[3], ".2d.syn_fits.txt"), header = F, sep = "\t") 
}
if (args[4] == "complex") 
{ 
	data = read.delim(paste0("../results/dadi_outputs/sim_", args[1], "-", args[2], "-", args[3], ".2d.complex_dem_syn_fits.txt"), header = F, sep = "\t") 
}

if (args[4] == "flexible") 
{ 
	data = read.delim(paste0("../results/dadi_outputs/sim_", args[1], "-", args[2], "-", args[3], ".2d.flexible_dem_syn_fits.txt"), header = F, sep = "\t") 
}

# print(data)
data = subset(data, V1 != "Inf")
best_fit = tail(data[order(data$V1),], n=1)
best_fit = cbind(best_fit, 2)
print(best_fit)


if (args[4] == "true")    { write.table(best_fit, file = paste0("../results/dadi_outputs/sim_", args[1], "-", args[2], "-", args[3], ".2d.syn_best_fit.txt"), sep = "\t", quote = F, col.names = F, row.names = F, append = F) }
if (args[4] == "complex") { write.table(best_fit, file = paste0("../results/dadi_outputs/sim_", args[1], "-", args[2], "-", args[3], ".2d.complex_dem_syn_best_fit.txt"), sep = "\t", quote = F, col.names = F, row.names = F, append = F) }
if (args[4] == "flexible") { write.table(best_fit, file = paste0("../results/dadi_outputs/sim_", args[1], "-", args[2], "-", args[3], ".2d.flexible_dem_syn_best_fit.txt"), sep = "\t", quote = F, col.names = F, row.names = F, append = F) }
