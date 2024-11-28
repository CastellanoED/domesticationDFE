source("postprocessing.R")
args = commandArgs(trailingOnly=TRUE)
library(dplyr)

args[1] #= "0" #mig
args[2] #= "2" #possel
args[3] #= "0" #change
args[4] #= "1" #iter

###########################
# parse the polyDFE output
###########################

# multiple runs of polyDFE can be stored in the same list
# for easier processing
est = list()
polyDFEfiles = c(paste0("../results/polyDFE_outputs/sim_", args[1], "-", args[2], "-", args[3], "-", args[4], "_1_polyDFE.txt"),paste0("../results/polyDFE_outputs/sim_", args[1], "-", args[2], "-", args[3], "-", args[4], "_10_polyDFE.txt"),paste0("../results/polyDFE_outputs/sim_", args[1], "-", args[2], "-", args[3], "-", args[4], "_2_polyDFE.txt"),paste0("../results/polyDFE_outputs/sim_", args[1], "-", args[2], "-", args[3], "-", args[4], "_20_polyDFE.txt"),paste0("../results/polyDFE_outputs/sim_", args[1], "-", args[2], "-", args[3], "-", args[4], "_30_polyDFE.txt"))

for (filename in polyDFEfiles)
{
    est = c(est, parseOutput(filename))
}
print(length(est))

# what are the gradients?
grad = sapply(est, function(e) e$criteria)
print(grad)
# for which files are the gradients a bit large?
polyDFEfiles[which(grad > 0.01)]
# for those runs it might help to run some basin hopping iterations on top

# on what input files was polyDFE ran?
print(sapply(est, function(e) e$input))

#print(estimateAlpha(est[[5]]))


###########################
## alpha from DFE w/o substitutions
###########################
alpha =  data.frame(t(sapply(est, function(e) c(estimateAlpha(e, supLimit = 1) ))))
print(alpha)
alpha$model = c(1, 10, 2, 20, 30)


###########################
## summarize the DFE
###########################
# we can also change the binning
dfe = data.frame(t(sapply(est, getDiscretizedDFE, sRanges = c(-10,-1, 0, 1, 10))))
dfe_wild = select(dfe, X1, X3, X5, X7, X9, X11)
dfe_dom  = select(dfe, X2, X4, X6, X8, X10, X12)

colnames(dfe_wild) <- c("< -10", "(-10, -1)", " (-1, 0)", "(0, 1)", "(1, 10)", "10 <")
colnames(dfe_dom) <- c("< -10", "(-10, -1)", " (-1, 0)", "(0, 1)", "(1, 10)", "10 <")
dfe_wild$model = c(1, 10, 2, 20, 30)
dfe_dom$model = c(1, 10, 2, 20, 30)

dfe_wild$population = "Wild" 
dfe_dom$population = "Domestic"


###########################
# model testing
###########################
# if only est1 is given, it just returns the AIC
aic = compareModels(est)
print(aic)
print(polyDFEfiles[which.min(aic$AIC[, "AIC"])])
#print(compareModels(est))
#1, 10, 2, 20
print(compareModels(est[1], est[2])) #1 vs 10 shared shape delDFE?
print(compareModels(est[1], est[3])) #1 vs 2 beneficial mutations with independent shape?
print(compareModels(est[3], est[4])) #2 vs 20 shared shape fullDFE?  
print(compareModels(est[2], est[4])) #10 vs 20 beneficial mutations with shared shape?
print(compareModels(est[4], est[5])) #20 vs 30 shared beneficial mutations with shared shape or independent beneficial mutations with shared shape?

LRT1 = data.frame(compareModels(est[1], est[2]))
LRT2 = data.frame(compareModels(est[1], est[3]))
LRT3 = data.frame(compareModels(est[3], est[4]))
LRT4 = data.frame(compareModels(est[2], est[4]))
LRT5 = data.frame(compareModels(est[4], est[5]))

LRT = cbind(LRT1$LRT.p.value, LRT2$LRT.p.value, LRT3$LRT.p.value, LRT4$LRT.p.value, LRT5$LRT.p.value, args[1], args[2], args[3], args[4])
write.table(LRT,file="../results/data_frames/polyDFE_lrt_replicate.txt", row.names = F, col.names = F, quote = F, append = TRUE)


###########################
# model averaging
###########################
# calculate the AIC weights
aic_weights = getAICweights(est)
print(aic_weights)

# calculate model-averaged DFE
estimates1 = data.frame(t(data.frame(est[[1]]$values)))
estimates2 = data.frame(t(data.frame(est[[2]]$values)))
estimates3 = data.frame(t(data.frame(est[[3]]$values)))
estimates4 = data.frame(t(data.frame(est[[4]]$values)))
estimates5 = data.frame(t(data.frame(est[[5]]$values)))

S_d   = rbind(estimates1$S_d,estimates2$S_d,estimates3$S_d,estimates4$S_d,estimates5$S_d)
shape = rbind(estimates1$b,estimates2$b,estimates3$b,estimates4$b,estimates5$b)
S_b   = rbind(estimates1$S_b,estimates2$S_b,estimates3$S_b,estimates4$S_b,estimates5$S_b)
p_b   = rbind(estimates1$p_b,estimates2$p_b,estimates3$p_b,estimates4$p_b,estimates5$p_b)
theta   = rbind(estimates1$theta_bar,estimates2$theta_bar,estimates3$theta_bar,estimates4$theta_bar,estimates5$theta_bar)

aic_weights = data.frame(aic_weights)

### PRINT OUTPUTS ###

output = cbind(S_d,shape,S_b,p_b,theta,aic_weights,args[1],args[2], args[3], args[4])  
output$model = c(1, 10, 2, 20, 30)
write.table(output,file="../results/data_frames/polyDFE_parameters_replicates_allmodels.txt",row.names = F, col.names = F, quote = F, append = TRUE)

dfe_wild = cbind(dfe_wild, aic_weights, args[1], args[2], args[3], args[4])
dfe_dom = cbind(dfe_dom, aic_weights, args[1], args[2], args[3], args[4])
dfe_both = rbind(dfe_wild, dfe_dom)
write.table(dfe_both,file="../results/data_frames/polyDFE_DFE_replicates_allmodels.txt",row.names = F, col.names = F, quote = F, append = TRUE)

alpha = cbind(alpha, aic_weights, args[1], args[2], args[3], args[4])
write.table(alpha,file="../results/data_frames/polyDFE_alpha_replicates_allmodels.txt",row.names = F, col.names = F, quote = F, append = TRUE)
