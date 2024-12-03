library(dplyr)
library(ggplot2)
# library(ggpubr)

### LRT ###
data = read.delim("../results/data_frames/polyDFE_lrt_replicate.txt", sep = " ", header = F)
data = unique(data)
colnames(data) <- c("pvalue_model1_vs_model10", "pvalue_model1_vs_model2", "pvalue_model2_vs_model20", "pvalue_model10_vs_model20", "pvalue_model20_vs_model30", "mig", "possel", "pchange", "it")
data$params = paste0(data$mig, "-", data$possel, "-", data$pchange)
table(data$params) #number of replicates per scenario

data$params <- factor(data$params, levels = c('0-0-0','0-0-0.05', '0-0-0.25','0.01-0-0','0.01-0-0.05', '0.01-0-0.25', '0-2-0','0-2-0.05', '0-2-0.25', '0.01-2-0','0.01-2-0.05', '0.01-2-0.25', '0-20-0','0-20-0.05', '0-20-0.25', '0.01-20-0','0.01-20-0.05', '0.01-20-0.25', '0-200-0','0-200-0.05', '0-200-0.25', '0.01-200-0','0.01-200-0.05', '0.01-200-0.25'), labels = c("Mig=0%, Sb=0, pb=0%, Change=0%", "Mig=0%, Sb=0, pb=0%, Change=5%", "Mig=0%, Sb=0, pb=0%, Change=25%", "Mig=1%, Sb=0, pb=0%, Change=0%", "Mig=1%, Sb=0, pb=0%, Change=5%", "Mig=1%, Sb=0, pb=0%, Change=25%", "Mig=0%, Sb=1, pb=10%, Change=0%", "Mig=0%, Sb=1, pb=10%, Change=5%", "Mig=0%, Sb=1, pb=10%, Change=25%", "Mig=1%, Sb=1, pb=10%, Change=0%", "Mig=1%, Sb=1, pb=10%, Change=5%", "Mig=1%, Sb=1, pb=10%, Change=25%", "Mig=0%, Sb=10, pb=1%, Change=0%", "Mig=0%, Sb=10, pb=1%, Change=5%", "Mig=0%, Sb=10, pb=1%, Change=25%", "Mig=1%, Sb=10, pb=1%, Change=0%", "Mig=1%, Sb=10, pb=1%, Change=5%", "Mig=1%, Sb=10, pb=1%, Change=25%", "Mig=0%, Sb=100, pb=0.1%, Change=0%", "Mig=0%, Sb=100, pb=0.1%, Change=5%", "Mig=0%, Sb=100, pb=0.1%, Change=25%", "Mig=1%, Sb=100, pb=0.1%, Change=0%", "Mig=1%, Sb=100, pb=0.1%, Change=5%", "Mig=1%, Sb=100, pb=0.1%, Change=25%"), ordered = TRUE)

mean_pvalue = data.frame(rbind(tapply(data$pvalue_model1_vs_model10, data$params, mean),
tapply(data$pvalue_model1_vs_model2, data$params, mean),
tapply(data$pvalue_model2_vs_model20, data$params, mean),
tapply(data$pvalue_model10_vs_model20, data$params, mean),
tapply(data$pvalue_model20_vs_model30, data$params, mean)))
write.table(mean_pvalue, file="../results/data_frames/polyDFE_mean_pvalue_LRT.txt", quote = F, sep = "\t", append = FALSE)

median_pvalue = data.frame(rbind(tapply(data$pvalue_model1_vs_model10, data$params, median),
tapply(data$pvalue_model1_vs_model2, data$params, median),
tapply(data$pvalue_model2_vs_model20, data$params, median),
tapply(data$pvalue_model10_vs_model20, data$params, median),
tapply(data$pvalue_model20_vs_model30, data$params, median)))
write.table(median_pvalue, file="../results/data_frames/polyDFE_median_pvalue_LRT.txt", quote = F, sep = "\t", append = FALSE)


### DFE parameters ###
library(tidyr)
library(ggplot2)

data = read.delim("../results/data_frames/polyDFE_parameters_replicates_allmodels.txt", sep = " ", header = F)
colnames(data) <- c("Wild_Sd", "Domesticated_Sd", "Wild_b", "Domesticated_b", "Wild_Sb", "Domesticated_Sb", "Wild_pb", "Domesticated_pb", "Wild_theta", "Domesticated_theta", "df", "loglk", "deltaAIC", "AIC_weight", "mig", "possel", "pchange", "it", "Model")
data$params = paste0(data$mig, "-", data$possel, "-", data$pchange)

mutation_rate = 1.25e-7 * 2

data$Domesticated = data$Domesticated_theta / (4*mutation_rate)
data$Wild = data$Wild_theta / (4*mutation_rate)
data_mini = dplyr::select(data, "params", "Wild","Domesticated")
demog_Ne <- gather(data_mini, Population, Ne, Wild:Domesticated, factor_key=TRUE)
table(demog_Ne$params)

demog_Ne = data.frame(tapply(demog_Ne$Ne, list(demog_Ne$params, demog_Ne$Population), mean))
demog_Ne$params = row.names(demog_Ne)
demog_Ne <- gather(demog_Ne, Population, Ne, Wild:Domesticated, factor_key=TRUE)
demog_Ne$Method = "polyDFE"
write.table(demog_Ne, file="../results/data_frames/polyDFE_Ne.txt", quote = F, sep = "\t", row.names = F, append = FALSE)

scaling_domesticated = data.frame(tapply(data$Domesticated, data$params, mean))
scaling_wild = data.frame(tapply(data$Wild, data$params, mean))
scaling = cbind(scaling_wild, scaling_domesticated)
scaling$params = row.names(scaling)
colnames(scaling) <- c("Wild", "Domesticated", "params")
write.table(scaling, file="../results/data_frames/polyDFE_diversity_Ne_scaling.txt", quote = F, sep = "\t", row.names = F, append = FALSE)

data$w_Domesticated_Sd = data$AIC_weight * (data$Domesticated_Sd / data$Domesticated_theta * mutation_rate)
data$w_Domesticated_b = data$AIC_weight * data$Domesticated_b
data$w_Domesticated_Sb = data$AIC_weight * (data$Domesticated_Sb / data$Domesticated_theta * mutation_rate) 
data$w_Domesticated_pb = data$AIC_weight * data$Domesticated_pb 

Domesticated_Sd = data.frame(tapply(data$w_Domesticated_Sd, list(data$it, data$params), sum) ) 
Domesticated_b = data.frame(tapply(data$w_Domesticated_b, list(data$it, data$params), sum) ) 
Domesticated_Sb = data.frame(tapply(data$w_Domesticated_Sb, list(data$it, data$params), sum) ) 
Domesticated_pb = data.frame(tapply(data$w_Domesticated_pb, list(data$it, data$params), sum) ) 

data$w_Wild_Sd = data$AIC_weight * (data$Wild_Sd / data$Wild_theta * mutation_rate)
data$w_Wild_b = data$AIC_weight * data$Wild_b
data$w_Wild_Sb = data$AIC_weight * (data$Wild_Sb / data$Wild_theta * mutation_rate) 
data$w_Wild_pb = data$AIC_weight * data$Wild_pb 

Wild_Sd = data.frame(tapply(data$w_Wild_Sd, list(data$it, data$params), sum) ) 
Wild_b = data.frame(tapply(data$w_Wild_b, list(data$it, data$params), sum) ) 
Wild_Sb = data.frame(tapply(data$w_Wild_Sb, list(data$it, data$params), sum) ) 
Wild_pb = data.frame(tapply(data$w_Wild_pb, list(data$it, data$params), sum) ) 

data_long_Domesticated_Sd <- gather(Domesticated_Sd, params, factor_key=TRUE)
data_long_Domesticated_b <- gather(Domesticated_b, params, factor_key=TRUE)
data_long_Domesticated_Sb <- gather(Domesticated_Sb, params, factor_key=TRUE)
data_long_Domesticated_pb <- gather(Domesticated_pb, params, factor_key=TRUE)

data_long_Wild_Sd <- gather(Wild_Sd, params, factor_key=TRUE)
data_long_Wild_b <- gather(Wild_b, params, factor_key=TRUE)
data_long_Wild_Sb <- gather(Wild_Sb, params, factor_key=TRUE)
data_long_Wild_pb <- gather(Wild_pb, params, factor_key=TRUE)

colnames(data_long_Domesticated_Sd) <- c("Simulation", "Domesticated")
colnames(data_long_Domesticated_b) <- c("Simulation", "Domesticated")
colnames(data_long_Domesticated_Sb) <- c("Simulation", "Domesticated")
colnames(data_long_Domesticated_pb) <- c("Simulation", "Domesticated")

colnames(data_long_Wild_Sd) <- c("Simulation", "Wild")
colnames(data_long_Wild_b) <- c("Simulation", "Wild")
colnames(data_long_Wild_Sb) <- c("Simulation", "Wild")
colnames(data_long_Wild_pb) <- c("Simulation", "Wild")

data_long_Domesticated_Sd$Replicate <- rownames(data_long_Domesticated_Sd)
data_long_Domesticated_b$Replicate <- rownames(data_long_Domesticated_b)
data_long_Domesticated_Sb$Replicate <- rownames(data_long_Domesticated_Sb)
data_long_Domesticated_pb$Replicate <- rownames(data_long_Domesticated_pb)

data_long_Wild_Sd$Replicate <- rownames(data_long_Wild_Sd)
data_long_Wild_b$Replicate <- rownames(data_long_Wild_b)
data_long_Wild_Sb$Replicate <- rownames(data_long_Wild_Sb)
data_long_Wild_pb$Replicate <- rownames(data_long_Wild_pb)

data_wide_Sd = merge(data_long_Domesticated_Sd, data_long_Wild_Sd)
data_long_Sd <- gather(data_wide_Sd, Population, Sd, Domesticated:Wild, factor_key=TRUE)
summary(data_long_Sd)

data_long_stat = data_long_Sd
data_long_stat$stat = data_long_stat$Sd
  
data_long_stat = na.omit(data_long_stat)

median_stat = data.frame(t(round(tapply(data_long_stat$stat, list(data_long_stat$Population, data_long_stat$Simulation), quantile, p=c(0.5) ), digits=2)))
ci5_stat = data.frame(t(round(tapply(data_long_stat$stat, list(data_long_stat$Population, data_long_stat$Simulation), quantile, p=c(0.05) ), digits=2)))
ci95_stat = data.frame(t(round(tapply(data_long_stat$stat, list(data_long_stat$Population, data_long_stat$Simulation), quantile, p=c(0.95) ), digits=2)))

wild_col = paste0(median_stat$Wild, "(", ci5_stat$Wild, ", ", ci95_stat$Wild, ")")
domesticated_col = paste0(median_stat$Domesticated, "(", ci5_stat$Domesticated, ", ", ci95_stat$Domesticated, ")")
sim_col = row.names(median_stat)
stats1 = data.frame(cbind(sim_col, wild_col, domesticated_col))
write.table(stats1, file="../results/data_frames/polyDFE_AIC_w_Sd.txt", quote = F, sep = "\t", row.names = F)
  

data_wide_b = merge(data_long_Domesticated_b, data_long_Wild_b)
data_long_b <- gather(data_wide_b, Population, b, Domesticated:Wild, factor_key=TRUE)
summary(data_long_b)

data_long_stat = data_long_b
data_long_stat$stat = data_long_stat$b

data_long_stat = na.omit(data_long_stat)

median_stat = data.frame(t(round(tapply(data_long_stat$stat, list(data_long_stat$Population, data_long_stat$Simulation), quantile, p=c(0.5) ), digits=2)))
ci5_stat = data.frame(t(round(tapply(data_long_stat$stat, list(data_long_stat$Population, data_long_stat$Simulation), quantile, p=c(0.05) ), digits=2)))
ci95_stat = data.frame(t(round(tapply(data_long_stat$stat, list(data_long_stat$Population, data_long_stat$Simulation), quantile, p=c(0.95) ), digits=2)))

wild_col = paste0(median_stat$Wild, "(", ci5_stat$Wild, ", ", ci95_stat$Wild, ")")
domesticated_col = paste0(median_stat$Domesticated, "(", ci5_stat$Domesticated, ", ", ci95_stat$Domesticated, ")")
sim_col = row.names(median_stat)
stats1 = data.frame(cbind(sim_col, wild_col, domesticated_col))
write.table(stats1, file="../results/data_frames/polyDFE_AIC_w_b.txt", quote = F, sep = "\t", row.names = F)


data_wide_Sb = merge(data_long_Domesticated_Sb, data_long_Wild_Sb)
data_long_Sb <- gather(data_wide_Sb, Population, Sb, Domesticated:Wild, factor_key=TRUE)
summary(data_long_Sb)

data_long_stat = data_long_Sb
data_long_stat$stat = data_long_stat$Sb

data_long_stat = na.omit(data_long_stat)

median_stat = data.frame(t(round(tapply(data_long_stat$stat, list(data_long_stat$Population, data_long_stat$Simulation), quantile, p=c(0.5) ), digits=2)))
ci5_stat = data.frame(t(round(tapply(data_long_stat$stat, list(data_long_stat$Population, data_long_stat$Simulation), quantile, p=c(0.05) ), digits=2)))
ci95_stat = data.frame(t(round(tapply(data_long_stat$stat, list(data_long_stat$Population, data_long_stat$Simulation), quantile, p=c(0.95) ), digits=2)))

wild_col = paste0(median_stat$Wild, "(", ci5_stat$Wild, ", ", ci95_stat$Wild, ")")
domesticated_col = paste0(median_stat$Domesticated, "(", ci5_stat$Domesticated, ", ", ci95_stat$Domesticated, ")")
sim_col = row.names(median_stat)
stats1 = data.frame(cbind(sim_col, wild_col, domesticated_col))
write.table(stats1, file="../results/data_frames/polyDFE_AIC_w_Sb.txt", quote = F, sep = "\t", row.names = F)


data_wide_pb = merge(data_long_Domesticated_pb, data_long_Wild_pb)
data_long_pb <- gather(data_wide_pb, Population, pb, Domesticated:Wild, factor_key=TRUE)
summary(data_long_pb)

data_long_stat = data_long_pb
data_long_stat$stat = data_long_stat$pb

data_long_stat = na.omit(data_long_stat)

median_stat = data.frame(t(round(tapply(data_long_stat$stat, list(data_long_stat$Population, data_long_stat$Simulation), quantile, p=c(0.5) ), digits=4)))
ci5_stat = data.frame(t(round(tapply(data_long_stat$stat, list(data_long_stat$Population, data_long_stat$Simulation), quantile, p=c(0.05) ), digits=4)))
ci95_stat = data.frame(t(round(tapply(data_long_stat$stat, list(data_long_stat$Population, data_long_stat$Simulation), quantile, p=c(0.95) ), digits=4)))

wild_col = paste0(median_stat$Wild, "(", ci5_stat$Wild, ", ", ci95_stat$Wild, ")")
domesticated_col = paste0(median_stat$Domesticated, "(", ci5_stat$Domesticated, ", ", ci95_stat$Domesticated, ")")
sim_col = row.names(median_stat)
stats1 = data.frame(cbind(sim_col, wild_col, domesticated_col))
write.table(stats1, file="../results/data_frames/polyDFE_AIC_w_pb.txt", quote = F, sep = "\t", row.names = F)


polyDFE_AIC_weighted_params1 = inner_join(data_long_Sd, data_long_b)
polyDFE_AIC_weighted_params2 = inner_join(data_long_Sb, data_long_pb)
polyDFE_AIC_weighted_params = inner_join(polyDFE_AIC_weighted_params1, polyDFE_AIC_weighted_params2)
write.table(polyDFE_AIC_weighted_params, file = "../results/data_frames/polyDFE_best_fits.txt", row.names = F, col.names = T, sep = "\t", quote = F)


