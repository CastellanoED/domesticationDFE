library(dplyr)
library(ggplot2)
# library(ggpubr)

setwd("~/Desktop/xdisk_dcastellano/domestication_DFE/scripts")

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


myplot = ggplot(data, aes(x=log10(pvalue_model1_vs_model10), fill=as.factor(params))) +
  ggtitle("LRT: delDFE with independent shapes vs\n delDFE with shared shapes") +
  #geom_density(adjust=2, alpha=0.5) +
  geom_boxplot() +
  facet_wrap(~ as.factor(params), ncol = 3) +
  scale_fill_manual(values=c("grey", "grey", "grey","grey", "grey", "grey","grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey")) +
  #scale_colour_manual(values=cbbPalette) +
  geom_vline(aes(xintercept=log10(0.05)),
             color="black", linetype="dotted", size=1) +
  #xlim(0,1) +
  theme_classic() +
  #theme(panel.spacing.x = unit(0, "mm")) +
  xlab("log10(p-value)") +
  theme(text = element_text(size=15), legend.position="none",axis.text.y=element_blank(),axis.title.y=element_blank())
ggsave(myplot, file= "../plots/polyDFE/LRT1.png", width=12, height=15.7, dpi=100)


myplot = ggplot(data, aes(x=log10(pvalue_model1_vs_model2), fill=as.factor(params))) +
  ggtitle("LRT: fullDFE with independent shapes vs\n delDFE with independent shapes") +
  #geom_density(adjust=2, alpha=0.5) +
  geom_boxplot() +
  facet_wrap(~ as.factor(params), ncol = 3) +
  scale_fill_manual(values=c("grey", "grey", "grey","grey", "grey", "grey","grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey")) +
  #scale_colour_manual(values=cbbPalette) +
  geom_vline(aes(xintercept=log10(0.05)),
             color="black", linetype="dotted", size=1) +
  #xlim(0,1) +
  theme_classic() +
  #theme(panel.spacing.x = unit(0, "mm")) +
  xlab("log10(p-value)") +
  theme(text = element_text(size=15), legend.position="none",axis.text.y=element_blank(),axis.title.y=element_blank())
ggsave(myplot, file= "../plots/polyDFE/LRT2.png", width=12, height=15.7, dpi=100)


myplot = ggplot(data, aes(x=log10(pvalue_model2_vs_model20), fill=as.factor(params))) +
  ggtitle("LRT: fullDFE with independent shapes vs\n fullDFE with shared shapes") +
  #geom_density(adjust=2, alpha=0.5) +
  geom_boxplot() +
  facet_wrap(~ as.factor(params), ncol = 3) +
  scale_fill_manual(values=c("grey", "grey", "grey","grey", "grey", "grey","grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey")) +
  #scale_colour_manual(values=cbbPalette) +
  geom_vline(aes(xintercept=log10(0.05)),
             color="black", linetype="dotted", size=1) +
  #xlim(0,1) +
  theme_classic() +
  #theme(panel.spacing.x = unit(0, "mm")) +
  xlab("log10(p-value)") +
  theme(text = element_text(size=15), legend.position="none",axis.text.y=element_blank(),axis.title.y=element_blank())
ggsave(myplot, file= "../plots/polyDFE/LRT3.png", width=12, height=15.7, dpi=100)


myplot = ggplot(data, aes(x=log10(pvalue_model10_vs_model20), fill=as.factor(params) ) ) +
  ggtitle("LRT: delDFE with shared shapes vs\n fullDFE with shared shapes") +
  geom_boxplot() +
  facet_wrap(~ as.factor(params), ncol = 3) +
  scale_fill_manual(values=c("grey", "grey", "grey","grey", "grey", "grey","grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey")) +
  #scale_colour_manual(values=cbbPalette) +
  geom_vline(aes(xintercept=log10(0.05)),
             color="black", linetype="dotted", size=1) +
  #xlim(0,1) +
  theme_classic() +
  #theme(panel.spacing.x = unit(0, "mm")) +
  xlab("log10(p-value)") +
  theme(text = element_text(size=15), legend.position="none",axis.text.y=element_blank(),axis.title.y=element_blank())
ggsave(myplot, file= "../plots/polyDFE/LRT4.png", width=12, height=15.7, dpi=100)


myplot = ggplot(data, aes(x=log10(pvalue_model20_vs_model30), fill=as.factor(params))) +
  ggtitle("LRT: fullDFE with shared shapes vs \nfullDFE with shared shapes and shared positive DFE") +
  geom_boxplot() +
  facet_wrap(~ as.factor(params), ncol = 3) +
  scale_fill_manual(values=c("grey", "grey", "grey","grey", "grey", "grey","grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey","grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey")) +
  #scale_colour_manual(values=cbbPalette) +
  geom_vline(aes(xintercept=log10(0.05)),
             color="black", linetype="dotted", size=1) +
  #xlim(0,1) +
  theme_classic() +
  #theme(panel.spacing.x = unit(0, "mm")) +
  xlab("log10(p-value)") +
  theme(text = element_text(size=15), legend.position="none",axis.text.y=element_blank(),axis.title.y=element_blank())
ggsave(myplot, file= "../plots/polyDFE/LRT5.png", width=12, height=15.7, dpi=100)




### DFE parameters ###
library(tidyr)
library(ggplot2)

data = read.delim("../results/data_frames/polyDFE_parameters_replicates_allmodels.txt", sep = " ", header = F)
colnames(data) <- c("Wild_Sd", "Domesticated_Sd", "Wild_b", "Domesticated_b", "Wild_Sb", "Domesticated_Sb", "Wild_pb", "Domesticated_pb", "Wild_theta", "Domesticated_theta", "df", "loglk", "deltaAIC", "AIC_weight", "mig", "possel", "pchange", "it", "Model")
data$params = paste0(data$mig, "-", data$possel, "-", data$pchange)

mutation_rate = 1.25e-7 * 2

4*5000*mutation_rate
4*10000*mutation_rate/2

  data$Domesticated = data$Domesticated_theta / (4*mutation_rate)
  data$Wild = data$Wild_theta / (4*mutation_rate)
  data_mini = dplyr::select(data, "params", "Wild","Domesticated")
  demog_Ne <- gather(data_mini, Population, Ne, Wild:Domesticated, factor_key=TRUE)
  table(demog_Ne$params)
  #demog_Ne$params <- factor(demog_Ne$params, levels = c('0-0-0','0-0-0.05', '0-0-0.25','0.01-0-0','0.01-0-0.05', '0.01-0-0.25', '0-2-0','0-2-0.05', '0-2-0.25', '0.01-2-0','0.01-2-0.05', '0.01-2-0.25', '0-20-0','0-20-0.05', '0-20-0.25', '0.01-20-0','0.01-20-0.05', '0.01-20-0.25', '0-200-0','0-200-0.05', '0-200-0.25', '0.01-200-0','0.01-200-0.05', '0.01-200-0.25'), labels = c("Mig=0%, Sb=0, pb=0%, Change=0%", "Mig=0%, Sb=0, pb=0%, Change=5%", "Mig=0%, Sb=0, pb=0%, Change=25%", "Mig=1%, Sb=0, pb=0%, Change=0%", "Mig=1%, Sb=0, pb=0%, Change=5%", "Mig=1%, Sb=0, pb=0%, Change=25%", "Mig=0%, Sb=1, pb=10%, Change=0%", "Mig=0%, Sb=1, pb=10%, Change=5%", "Mig=0%, Sb=1, pb=10%, Change=25%", "Mig=1%, Sb=1, pb=10%, Change=0%", "Mig=1%, Sb=1, pb=10%, Change=5%", "Mig=1%, Sb=1, pb=10%, Change=25%", "Mig=0%, Sb=10, pb=1%, Change=0%", "Mig=0%, Sb=10, pb=1%, Change=5%", "Mig=0%, Sb=10, pb=1%, Change=25%", "Mig=1%, Sb=10, pb=1%, Change=0%", "Mig=1%, Sb=10, pb=1%, Change=5%", "Mig=1%, Sb=10, pb=1%, Change=25%", "Mig=0%, Sb=100, pb=0.1%, Change=0%", "Mig=0%, Sb=100, pb=0.1%, Change=5%", "Mig=0%, Sb=100, pb=0.1%, Change=25%", "Mig=1%, Sb=100, pb=0.1%, Change=0%", "Mig=1%, Sb=100, pb=0.1%, Change=5%", "Mig=1%, Sb=100, pb=0.1%, Change=25%"), ordered = TRUE)

  myplot = ggplot(demog_Ne, aes(x=Ne, fill=Population)) +
    ggtitle("Ne best demographic model vs polyDFE") +
    geom_density() +
    facet_wrap(~ params, ncol = 3) +
    #scale_fill_grey() +
    geom_vline(xintercept = 5000, linetype="dashed") +
    xlim(0,10000) +
    theme_classic() +
    xlab("Ne") +
    theme(text = element_text(size=15), axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position="bottom")
  ggsave(myplot, file="../plots/polyDFE/Ne.png", width=12, height=15.7, dpi=100)
  
  demog_Ne = data.frame(tapply(demog_Ne$Ne, list(demog_Ne$params, demog_Ne$Population), mean))
  demog_Ne$params = row.names(demog_Ne)
  demog_Ne <- gather(demog_Ne, Population, Ne, Wild:Domesticated, factor_key=TRUE)
  
  myplot = ggplot(demog_Ne, aes(x=Ne, y = 1, color=Population)) +
    ggtitle("Ne dadi's best demographic model vs polyDFE") +
    geom_point(size=5,alpha=0.75) +
    scale_color_brewer(palette="Dark2") +
    facet_wrap(~ params, ncol = 3) +
    #scale_fill_grey() +
    geom_vline(xintercept = 5000, linetype="dashed") +
    xlim(0,10000) +
    theme_classic() +
    xlab("Ne") +
    theme(text = element_text(size=15), axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position="bottom")
  ggsave(myplot, file="../plots/polyDFE/Ne.png", width=12, height=15.7, dpi=100)
  
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
  
##data_long_stat$Simulation <- factor(data_long_stat$Simulation, levels = c('X0.2.0', 'X0.2.0.05', 'X0.2.0.25', 'X0.01.2.0', 'X0.01.2.0.05', 'X0.01.2.0.25', 'X0.20.0', 'X0.20.0.05', 'X0.20.0.25', 'X0.01.20.0', 'X0.01.20.0.05', 'X0.01.20.0.25', 'X0.200.0', 'X0.200.0.05', 'X0.200.0.25', 'X0.01.200.0', 'X0.01.200.0.05', 'X0.01.200.0.25'), labels = c("ID = 1", "ID = 2", "ID = 3", "ID = 4", "ID = 5", "ID = 6", "ID = 7", "ID = 8", "ID = 9", "ID = 10", "ID = 11", "ID = 12", "ID = 13", "ID = 14", "ID = 15", "ID = 16", "ID = 17", "ID = 18"), ordered = TRUE)
data_long_stat$Simulation <- factor(data_long_stat$Simulation, levels = c('X0.0.0', 'X0.0.0.05', 'X0.0.0.25', 'X0.01.0.0', 'X0.01.0.0.05', 'X0.01.0.0.25',  'X0.2.0', 'X0.2.0.05', 'X0.2.0.25', 'X0.01.2.0', 'X0.01.2.0.05', 'X0.01.2.0.25', 'X0.20.0', 'X0.20.0.05', 'X0.20.0.25', 'X0.01.20.0', 'X0.01.20.0.05', 'X0.01.20.0.25', 'X0.200.0', 'X0.200.0.05', 'X0.200.0.25', 'X0.01.200.0', 'X0.01.200.0.05', 'X0.01.200.0.25'), labels = c("Mig=0%, Sb=0,\npb=0%, Change=0%", "Mig=0%, Sb=0,\npb=0%, Change=5%", "Mig=0%, Sb=0,\npb=0%, Change=25%", "Mig=1%, Sb=0,\npb=0%, Change=0%", "Mig=1%, Sb=0,\npb=0%, Change=5%", "Mig=1%, Sb=0,\npb=0%, Change=25%",    "Mig=0%, Sb=1,\npb=10%, Change=0%", "Mig=0%, Sb=1,\npb=10%, Change=5%", "Mig=0%, Sb=1,\npb=10%, Change=25%", "Mig=1%, Sb=1,\npb=10%, Change=0%", "Mig=1%, Sb=1,\npb=10%, Change=5%", "Mig=1%, Sb=1,\npb=10%, Change=25%", "Mig=0%, Sb=10,\npb=1%, Change=0%", "Mig=0%, Sb=10,\npb=1%, Change=5%", "Mig=0%, Sb=10,\npb=1%, Change=25%", "Mig=1%, Sb=10,\npb=1%, Change=0%", "Mig=1%, Sb=10,\npb=1%, Change=5%", "Mig=1%, Sb=10,\npb=1%, Change=25%", "Mig=0%, Sb=100,\npb=0.1%, Change=0%", "Mig=0%, Sb=100,\npb=0.1%, Change=5%", "Mig=0%, Sb=100,\npb=0.1%, Change=25%", "Mig=1%, Sb=100,\npb=0.1%, Change=0%", "Mig=1%, Sb=100,\npb=0.1%, Change=5%", "Mig=1%, Sb=100,\npb=0.1%, Change=25%"), ordered = TRUE)
##        demog_Ne$params <- factor(demog_Ne$params,           levels = c('0-2-0',   '0-2-0.05',  '0-2-0.25',  '0.01-2-0',  '0.01-2-0.05',  '0.01-2-0.25',  '0-20-0',  '0-20-0.05',  '0-20-0.25',  '0.01-20-0',  '0.01-20-0.05',  '0.01-20-0.25',  '0-200-0',  '0-200-0.05',  '0-200-0.25',  '0.01-200-0',  '0.01-200-0.05',  '0.01-200-0.25'), labels = c("Mig=0%, Sb=1, pb=10%, Change=0%", "Mig=0%, Sb=1, pb=10%, Change=5%", "Mig=0%, Sb=1, pb=10%, Change=25%", "Mig=1%, Sb=1, pb=10%, Change=0%", "Mig=1%, Sb=1, pb=10%, Change=5%", "Mig=1%, Sb=1, pb=10%, Change=25%", "Mig=0%, Sb=10, pb=1%, Change=0%", "Mig=0%, Sb=10, pb=1%, Change=5%", "Mig=0%, Sb=10, pb=1%, Change=25%", "Mig=1%, Sb=10, pb=1%, Change=0%", "Mig=1%, Sb=10, pb=1%, Change=5%", "Mig=1%, Sb=10, pb=1%, Change=25%", "Mig=0%, Sb=100, pb=0.1%, Change=0%", "Mig=0%, Sb=100, pb=0.1%, Change=5%", "Mig=0%, Sb=100, pb=0.1%, Change=25%", "Mig=1%, Sb=100, pb=0.1%, Change=0%", "Mig=1%, Sb=100, pb=0.1%, Change=5%", "Mig=1%, Sb=100, pb=0.1%, Change=25%"), ordered = TRUE)


myplot = ggplot(data_long_stat, aes(x=-stat, fill=Population)) +
  ggtitle("C") +
  geom_density(adjust=2, alpha=0.5) +
  #geom_boxplot() + 
  facet_wrap(~ Simulation, ncol = 3) +
  scale_fill_grey() +
  geom_vline(aes(xintercept=(100/10000)),
             color="black", linetype="dotted", size=1) +
  xlim(0,0.02) +
  theme_classic() +
  #theme(panel.spacing.x = unit(0, "mm")) +
  xlab("") +
  theme(text = element_text(size=25), axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position="none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
panelC=myplot
ggsave(myplot, file="../plots/polyDFE/Sd.png", width=12, height=15.7, dpi=100)



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

data_long_stat$Simulation <- factor(data_long_stat$Simulation, levels = c('X0.0.0', 'X0.0.0.05', 'X0.0.0.25', 'X0.01.0.0', 'X0.01.0.0.05', 'X0.01.0.0.25',  'X0.2.0', 'X0.2.0.05', 'X0.2.0.25', 'X0.01.2.0', 'X0.01.2.0.05', 'X0.01.2.0.25', 'X0.20.0', 'X0.20.0.05', 'X0.20.0.25', 'X0.01.20.0', 'X0.01.20.0.05', 'X0.01.20.0.25', 'X0.200.0', 'X0.200.0.05', 'X0.200.0.25', 'X0.01.200.0', 'X0.01.200.0.05', 'X0.01.200.0.25'), labels = c("Mig=0%, Sb=0,\npb=0%, Change=0%", "Mig=0%, Sb=0,\npb=0%, Change=5%", "Mig=0%, Sb=0,\npb=0%, Change=25%", "Mig=1%, Sb=0,\npb=0%, Change=0%", "Mig=1%, Sb=0,\npb=0%, Change=5%", "Mig=1%, Sb=0,\npb=0%, Change=25%",    "Mig=0%, Sb=1,\npb=10%, Change=0%", "Mig=0%, Sb=1,\npb=10%, Change=5%", "Mig=0%, Sb=1,\npb=10%, Change=25%", "Mig=1%, Sb=1,\npb=10%, Change=0%", "Mig=1%, Sb=1,\npb=10%, Change=5%", "Mig=1%, Sb=1,\npb=10%, Change=25%", "Mig=0%, Sb=10,\npb=1%, Change=0%", "Mig=0%, Sb=10,\npb=1%, Change=5%", "Mig=0%, Sb=10,\npb=1%, Change=25%", "Mig=1%, Sb=10,\npb=1%, Change=0%", "Mig=1%, Sb=10,\npb=1%, Change=5%", "Mig=1%, Sb=10,\npb=1%, Change=25%", "Mig=0%, Sb=100,\npb=0.1%, Change=0%", "Mig=0%, Sb=100,\npb=0.1%, Change=5%", "Mig=0%, Sb=100,\npb=0.1%, Change=25%", "Mig=1%, Sb=100,\npb=0.1%, Change=0%", "Mig=1%, Sb=100,\npb=0.1%, Change=5%", "Mig=1%, Sb=100,\npb=0.1%, Change=25%"), ordered = TRUE)

myplot = ggplot(data_long_stat, aes(x=(stat), fill=Population)) +
  ggtitle("A") +
  geom_density(adjust=2, alpha=0.5) +
  #geom_boxplot() + 
  facet_wrap(~ Simulation, ncol = 3) +
  scale_fill_grey() +
  geom_vline(aes(xintercept=(0.3)),
             color="black", linetype="dotted", size=1) +
  xlim(0,1.75) +
  theme_classic() +
  #theme(panel.spacing.x = unit(0, "mm")) +
  xlab("") +
  theme(text = element_text(size=25), axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position="none")
panelA = myplot
ggsave("../plots/polyDFE/b.png", width=12, height=15.7, dpi=100)


# library(gridExtra)
# myplot = grid.arrange(panelA, panelB, panelC, panelD, nrow = 2)
# ggsave(myplot, file="../plots/Figure3.png", width=20, height=24, dpi=500)




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

data_long_stat$Simulation <- factor(data_long_stat$Simulation, levels = c('X0.0.0', 'X0.0.0.05', 'X0.0.0.25', 'X0.01.0.0', 'X0.01.0.0.05', 'X0.01.0.0.25',  'X0.2.0', 'X0.2.0.05', 'X0.2.0.25', 'X0.01.2.0', 'X0.01.2.0.05', 'X0.01.2.0.25', 'X0.20.0', 'X0.20.0.05', 'X0.20.0.25', 'X0.01.20.0', 'X0.01.20.0.05', 'X0.01.20.0.25', 'X0.200.0', 'X0.200.0.05', 'X0.200.0.25', 'X0.01.200.0', 'X0.01.200.0.05', 'X0.01.200.0.25'), labels = c("Mig=0%, Sb=0,\npb=0%, Change=0%", "Mig=0%, Sb=0,\npb=0%, Change=5%", "Mig=0%, Sb=0,\npb=0%, Change=25%", "Mig=1%, Sb=0,\npb=0%, Change=0%", "Mig=1%, Sb=0,\npb=0%, Change=5%", "Mig=1%, Sb=0,\npb=0%, Change=25%",    "Mig=0%, Sb=1,\npb=10%, Change=0%", "Mig=0%, Sb=1,\npb=10%, Change=5%", "Mig=0%, Sb=1,\npb=10%, Change=25%", "Mig=1%, Sb=1,\npb=10%, Change=0%", "Mig=1%, Sb=1,\npb=10%, Change=5%", "Mig=1%, Sb=1,\npb=10%, Change=25%", "Mig=0%, Sb=10,\npb=1%, Change=0%", "Mig=0%, Sb=10,\npb=1%, Change=5%", "Mig=0%, Sb=10,\npb=1%, Change=25%", "Mig=1%, Sb=10,\npb=1%, Change=0%", "Mig=1%, Sb=10,\npb=1%, Change=5%", "Mig=1%, Sb=10,\npb=1%, Change=25%", "Mig=0%, Sb=100,\npb=0.1%, Change=0%", "Mig=0%, Sb=100,\npb=0.1%, Change=5%", "Mig=0%, Sb=100,\npb=0.1%, Change=25%", "Mig=1%, Sb=100,\npb=0.1%, Change=0%", "Mig=1%, Sb=100,\npb=0.1%, Change=5%", "Mig=1%, Sb=100,\npb=0.1%, Change=25%"), ordered = TRUE)

myplot = ggplot(data_long_stat, aes(x=log10(stat), fill=Population)) +
  ggtitle("AIC weighted sb bootstrap distributions") +
  geom_density(adjust=2, alpha=0.5) +
  #geom_boxplot() + 
  facet_wrap(~ Simulation, ncol = 3) +
  scale_fill_grey() +
  geom_vline(aes(xintercept=log10(1/10000)),
             color="black", linetype="dotted", size=1) +
  geom_vline(aes(xintercept=log10(10/10000)),
             color="black", linetype="dotted", size=1) +
  geom_vline(aes(xintercept=log10(100/10000)),
             color="black", linetype="dotted", size=1) +
  #xlim(0,100) +
  theme_classic() +
  #theme(panel.spacing.x = unit(0, "mm")) +
  xlab("log10(sb)") +
  theme(text = element_text(size=15), axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position = "bottom")
ggsave(plot = myplot, "../plots/polyDFE/sb.png", width=12, height=15.7, dpi=100)


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

data_long_stat$Simulation <- factor(data_long_stat$Simulation, levels = c('X0.0.0', 'X0.0.0.05', 'X0.0.0.25', 'X0.01.0.0', 'X0.01.0.0.05', 'X0.01.0.0.25',  'X0.2.0', 'X0.2.0.05', 'X0.2.0.25', 'X0.01.2.0', 'X0.01.2.0.05', 'X0.01.2.0.25', 'X0.20.0', 'X0.20.0.05', 'X0.20.0.25', 'X0.01.20.0', 'X0.01.20.0.05', 'X0.01.20.0.25', 'X0.200.0', 'X0.200.0.05', 'X0.200.0.25', 'X0.01.200.0', 'X0.01.200.0.05', 'X0.01.200.0.25'), labels = c("Mig=0%, Sb=0,\npb=0%, Change=0%", "Mig=0%, Sb=0,\npb=0%, Change=5%", "Mig=0%, Sb=0,\npb=0%, Change=25%", "Mig=1%, Sb=0,\npb=0%, Change=0%", "Mig=1%, Sb=0,\npb=0%, Change=5%", "Mig=1%, Sb=0,\npb=0%, Change=25%",    "Mig=0%, Sb=1,\npb=10%, Change=0%", "Mig=0%, Sb=1,\npb=10%, Change=5%", "Mig=0%, Sb=1,\npb=10%, Change=25%", "Mig=1%, Sb=1,\npb=10%, Change=0%", "Mig=1%, Sb=1,\npb=10%, Change=5%", "Mig=1%, Sb=1,\npb=10%, Change=25%", "Mig=0%, Sb=10,\npb=1%, Change=0%", "Mig=0%, Sb=10,\npb=1%, Change=5%", "Mig=0%, Sb=10,\npb=1%, Change=25%", "Mig=1%, Sb=10,\npb=1%, Change=0%", "Mig=1%, Sb=10,\npb=1%, Change=5%", "Mig=1%, Sb=10,\npb=1%, Change=25%", "Mig=0%, Sb=100,\npb=0.1%, Change=0%", "Mig=0%, Sb=100,\npb=0.1%, Change=5%", "Mig=0%, Sb=100,\npb=0.1%, Change=25%", "Mig=1%, Sb=100,\npb=0.1%, Change=0%", "Mig=1%, Sb=100,\npb=0.1%, Change=5%", "Mig=1%, Sb=100,\npb=0.1%, Change=25%"), ordered = TRUE)

myplot = ggplot(data_long_stat, aes(x=log10(stat), fill=Population)) +
  ggtitle("AIC weighted pb bootstrap distributions") +
  #geom_density(adjust=2, alpha=0.5) +
  geom_boxplot(alpha=0.5) + 
  facet_wrap(~ Simulation, ncol = 3) +
  scale_fill_grey() +
  #geom_vline(aes(xintercept=0),
  #           color="black", linetype="dotted", size=1) +
  geom_vline(aes(xintercept=log10(0.1)),
             color="black", linetype="dotted", size=1) +
  geom_vline(aes(xintercept=log10(0.01)),
             color="black", linetype="dotted", size=1) +
  geom_vline(aes(xintercept=log10(0.001)),
             color="black", linetype="dotted", size=1) +
  xlim(-3,0) +
  theme_classic() +
  #theme(panel.spacing.x = unit(0, "mm")) +
  xlab("log10(pb)") +
  theme(text = element_text(size=15), axis.text.y=element_blank(), axis.title.y=element_blank(), legend.position = "bottom")
ggsave(plot = myplot, "../plots/polyDFE/pb.png", width=12, height=15.7, dpi=100)


polyDFE_AIC_weighted_params1 = inner_join(data_long_Sd, data_long_b)
polyDFE_AIC_weighted_params2 = inner_join(data_long_Sb, data_long_pb)
polyDFE_AIC_weighted_params = inner_join(polyDFE_AIC_weighted_params1, polyDFE_AIC_weighted_params2)
write.table(polyDFE_AIC_weighted_params, file = "../results/data_frames/polyDFE_best_fits.txt", row.names = F, col.names = T, sep = "\t", quote = F)

### Merging data_long dataframes with AIC weighted parameters to infer customized s intervals ###

dt_polyDFE =  inner_join(data_long_Sd, data_long_b)
dt_polyDFE =  inner_join(dt_polyDFE, data_long_pb)

dt_int = dt_polyDFE
dt_int = dplyr::select(dt_int, -Replicate)

colnames(dt_int) <- c("params", "Population", "sd", "shape", "pb")
dt_int$sd = -dt_int$sd

#dt_int$sd = 2*20*((dt_int$shape*dt_int$scale)/dt_int$theta_nsyn)
dt_int$neu = pgamma(q=c(0.001), shape=dt_int$shape, scale=(dt_int$sd/dt_int$shape))
dt_int$strong_del = 1-pgamma(q=c(0.01), shape=dt_int$shape, scale=(dt_int$sd/dt_int$shape))
dt_int$weak_del = 1 - dt_int$neu - dt_int$strong_del
dt_int$ben = dt_int$pb

dt_int$BENEFICIAL = dt_int$ben
dt_int$NEUTRAL = (1-dt_int$ben)*dt_int$neu
dt_int$WEAKLY = (1-dt_int$ben)*dt_int$weak_del
dt_int$STRONG = (1-dt_int$ben)*dt_int$strong_del


dt_int = dplyr::select(dt_int, params, Population, BENEFICIAL:STRONG)
dt_int_long <- gather(dt_int, coeff, Value, BENEFICIAL:STRONG, factor_key=TRUE)

value =  tapply(dt_int_long$Value, paste0(dt_int_long$params,"_",dt_int_long$Population,"_", dt_int_long$coeff),  function(x) quantile(unlist(x), probs=c(0.5) ) ) 
CI5 =  tapply(dt_int_long$Value, paste0(dt_int_long$params,"_",dt_int_long$Population,"_", dt_int_long$coeff),  function(x) quantile(unlist(x), probs=c(0.025) ) ) 
CI95 =  tapply(dt_int_long$Value, paste0(dt_int_long$params,"_",dt_int_long$Population,"_", dt_int_long$coeff),  function(x) quantile(unlist(x), probs=c(0.975) ) ) 


dt_int = data.frame(cbind(value, CI5, CI95))
dt_int$info = row.names(dt_int)

dt_int = dt_int %>%
  separate(info, c("scenario", "Population", "coeff"), "_")

dt_int$coeff = gsub("BENEFICIAL", "> 0", dt_int$coeff)
dt_int$coeff = gsub("NEUTRAL", "(-0.001, 0]", dt_int$coeff)
dt_int$coeff = gsub("WEAKLY", "(-0.01, -0.001]", dt_int$coeff)
dt_int$coeff = gsub("STRONG", "< -0.01", dt_int$coeff)
dt_int$scenario <- factor(dt_int$scenario, levels = c('X0.0.0', 'X0.0.0.05', 'X0.0.0.25', 'X0.01.0.0', 'X0.01.0.0.05', 'X0.01.0.0.25',  'X0.2.0', 'X0.2.0.05', 'X0.2.0.25', 'X0.01.2.0', 'X0.01.2.0.05', 'X0.01.2.0.25', 'X0.20.0', 'X0.20.0.05', 'X0.20.0.25', 'X0.01.20.0', 'X0.01.20.0.05', 'X0.01.20.0.25', 'X0.200.0', 'X0.200.0.05', 'X0.200.0.25', 'X0.01.200.0', 'X0.01.200.0.05', 'X0.01.200.0.25'), labels = c("Mig=0%, Sb=0, pb=0%, Change=0%", "Mig=0%, Sb=0, pb=0%, Change=5%", "Mig=0%, Sb=0, pb=0%, Change=25%", "Mig=1%, Sb=0, pb=0%, Change=0%", "Mig=1%, Sb=0, pb=0%, Change=5%", "Mig=1%, Sb=0, pb=0%, Change=25%", "Mig=0%, Sb=1, pb=10%, Change=0%", "Mig=0%, Sb=1, pb=10%, Change=5%", "Mig=0%, Sb=1, pb=10%, Change=25%", "Mig=1%, Sb=1, pb=10%, Change=0%", "Mig=1%, Sb=1, pb=10%, Change=5%", "Mig=1%, Sb=1, pb=10%, Change=25%", "Mig=0%, Sb=10, pb=1%, Change=0%", "Mig=0%, Sb=10, pb=1%, Change=5%", "Mig=0%, Sb=10, pb=1%, Change=25%", "Mig=1%, Sb=10, pb=1%, Change=0%", "Mig=1%, Sb=10, pb=1%, Change=5%", "Mig=1%, Sb=10, pb=1%, Change=25%", "Mig=0%, Sb=100, pb=0.1%, Change=0%", "Mig=0%, Sb=100, pb=0.1%, Change=5%", "Mig=0%, Sb=100, pb=0.1%, Change=25%", "Mig=1%, Sb=100, pb=0.1%, Change=0%", "Mig=1%, Sb=100, pb=0.1%, Change=5%", "Mig=1%, Sb=100, pb=0.1%, Change=25%"), ordered = TRUE)

dt_int$scaling = 1
dt_int$Value = "Inferred"
dt_int$population = dt_int$Population
dt_int$Method = "polyDFE"

write.table(dt_int, file="../results/data_frames/discretized_DFE_polyDFE.txt", col.names = T, row.names = F, sep = "\t", quote = F, append = F)


### DFE intervals ###
library(tidyr)
library(ggplot2)
library(dplyr)
# library(ggpubr)
library(stringr)
# library(ggthemes)

dt <- read.delim("../results/data_frames/polyDFE_DFE_replicates_allmodels.txt", sep = " ", header = F)
colnames(dt) <- c("< -10", "(-10, -1)", "(-1, 0)", "(0, 1)", "(1, 10)", "10 <", "model", "population", "degrees_of_freedom", "loglk", "delta_AIC", "AIC_weight", "mig", "possel", "pchange", "replicate")
dt$scenario = paste0(dt$mig, "-", dt$possel, "-", dt$pchange)
dt$scenario <- factor(dt$scenario, levels = c('0-0-0', '0-0-0.05', '0-0-0.25', '0.01-0-0', '0.01-0-0.05', '0.01-0-0.25',   '0-2-0', '0-2-0.05', '0-2-0.25', '0.01-2-0', '0.01-2-0.05', '0.01-2-0.25', '0-20-0', '0-20-0.05', '0-20-0.25', '0.01-20-0', '0.01-20-0.05', '0.01-20-0.25', '0-200-0', '0-200-0.05', '0-200-0.25', '0.01-200-0', '0.01-200-0.05', '0.01-200-0.25'), labels = c(seq(24)), ordered = TRUE)
table(dt$scenario)

dt$neu = dt$`(-1, 0)` + dt$`(0, 1)`
dt$ben = dt$`(1, 10)` + dt$`10 <`

#dt_wild = subset(dt, population == "Wild")
#dt_wild$population = "Wild"
#dt_domestic = subset(dt, population == "Domestic")
#dt_domestic$population = "Domestic"

# dt_wild = subset(dt, population == "Domestic" & model == "20")
# dt_wild$population = "Wild"
# dt_domestic = subset(dt, population == "Wild" & model == "20")
# dt_domestic$population = "Domestic"

#dt = rbind(dt_wild, dt_domestic)

#dt$scenario = as.factor(dt$scenario)
#dt$scenario <- factor(df$scenario, levels = c("1_9", "2_8", "9_8", "1_1", "3_1", "6_4", "1_6", "5_3", "3_3", "8_6"))
#dt$scenario = as.character(dt$scenario)

dt$`< -10 w` = dt$`< -10` * dt$AIC_weight
dt$`(-10, -1) w` = dt$`(-10, -1)` * dt$AIC_weight
dt$`(-1, 0) w` = dt$`(-1, 0)` * dt$AIC_weight
dt$`(0, 1) w` = dt$`(0, 1)` * dt$AIC_weight
dt$`(1, 10) w` = dt$`(1, 10)` * dt$AIC_weight
dt$`10 < w` = dt$`10 <` * dt$AIC_weight
dt$`(-1, 0) w` = dt$neu * dt$AIC_weight
dt$`0 < w` = dt$ben * dt$AIC_weight

# dt$`< -10 w` = dt$`< -10` 
# dt$`(-10, -1) w` = dt$`(-10, -1)` 
# dt$`(-1, 0) w` = dt$`(-1, 0)` 
# dt$`(0, 1) w` = dt$`(0, 1)` 
# dt$`(1, 10) w` = dt$`(1, 10)` 
# dt$`10 < w` = dt$`10 <` 
# dt$`(-1, 0) w` = dt$neu 
# dt$`0 < w` = dt$ben 

selection_coef = data.frame(tapply(dt$`< -10 w`, list(dt$replicate, dt$scenario, dt$population), sum, na.rm=T))
selection_coef_long <- gather(selection_coef, scenario_population, value, factor_key=TRUE)
tapply(selection_coef_long$value, selection_coef_long$scenario_population, mean, na.rm=T)
selection_coef_long_split = data.frame(str_split_fixed(selection_coef_long$scenario_population, fixed("."), 2))
selection_coef_long$scenario = selection_coef_long_split$X1 
selection_coef_long$population = selection_coef_long_split$X2 
selection_coef_long$coeff = "strong_del"

strong_del = selection_coef_long

selection_coef = data.frame(tapply(dt$`(-10, -1) w`, list(dt$replicate, dt$scenario, dt$population), sum, na.rm=T))
selection_coef_long <- gather(selection_coef, scenario_population, value, factor_key=TRUE)
tapply(selection_coef_long$value, selection_coef_long$scenario_population, mean)
selection_coef_long_split = data.frame(str_split_fixed(selection_coef_long$scenario_population, fixed("."), 2))
selection_coef_long$scenario = selection_coef_long_split$X1 
selection_coef_long$population = selection_coef_long_split$X2 
selection_coef_long$coeff = "weakly_del"

weak_del = selection_coef_long

selection_coef = data.frame(tapply(dt$`(-1, 0) w`, list(dt$replicate, dt$scenario, dt$population), sum, na.rm=T))
selection_coef_long <- gather(selection_coef, scenario_population, value, factor_key=TRUE)
tapply(selection_coef_long$value, selection_coef_long$scenario_population, mean)
selection_coef_long_split = data.frame(str_split_fixed(selection_coef_long$scenario_population, fixed("."), 2))
selection_coef_long$scenario = selection_coef_long_split$X1 
selection_coef_long$population = selection_coef_long_split$X2 
selection_coef_long$coeff = "neu"

neu = selection_coef_long

selection_coef = data.frame(tapply(dt$`0 < w`, list(dt$replicate, dt$scenario, dt$population), sum, na.rm=T))
selection_coef_long <- gather(selection_coef, scenario_population, value, factor_key=TRUE)
tapply(selection_coef_long$value, selection_coef_long$scenario_population, mean)
selection_coef_long_split = data.frame(str_split_fixed(selection_coef_long$scenario_population, fixed("."), 2))
selection_coef_long$scenario = selection_coef_long_split$X1 
selection_coef_long$population = selection_coef_long_split$X2 
selection_coef_long$coeff = "ben"

ben = selection_coef_long


# # selection_coef = data.frame(tapply(dt$`(-1, 0) w`, list(dt$replicate, dt$scenario, dt$population), sum))
# # selection_coef_long <- gather(selection_coef, scenario_population, value, factor_key=TRUE)
# # tapply(selection_coef_long$value, selection_coef_long$scenario_population, mean)
# # selection_coef_long_split = data.frame(str_split_fixed(selection_coef_long$scenario_population, fixed("."), 2))
# # selection_coef_long$scenario = selection_coef_long_split$X1 
# # selection_coef_long$population = selection_coef_long_split$X2 
# # selection_coef_long$coeff = "(-1, 0)"

# # neu_del = selection_coef_long


# # selection_coef = data.frame(tapply(dt$`(0, 1) w`, list(dt$replicate, dt$scenario, dt$population), sum))
# # selection_coef_long <- gather(selection_coef, scenario_population, value, factor_key=TRUE)
# # tapply(selection_coef_long$value, selection_coef_long$scenario_population, mean)
# # selection_coef_long_split = data.frame(str_split_fixed(selection_coef_long$scenario_population, fixed("."), 2))
# # selection_coef_long$scenario = selection_coef_long_split$X1 
# # selection_coef_long$population = selection_coef_long_split$X2 
# # selection_coef_long$coeff = "(0, 1)"

# # neu_ben = selection_coef_long

# # selection_coef = data.frame(tapply(dt$`(1, 10) w`, list(dt$replicate, dt$scenario, dt$population), sum))
# # selection_coef_long <- gather(selection_coef, scenario_population, value, factor_key=TRUE)
# # tapply(selection_coef_long$value, selection_coef_long$scenario_population, mean)
# # selection_coef_long_split = data.frame(str_split_fixed(selection_coef_long$scenario_population, fixed("."), 2))
# # selection_coef_long$scenario = selection_coef_long_split$X1 
# # selection_coef_long$population = selection_coef_long_split$X2 
# # selection_coef_long$coeff = "(1, 10)"

# # weak_ben = selection_coef_long

# # selection_coef = data.frame(tapply(dt$`10 < w`, list(dt$replicate, dt$scenario, dt$population), sum))
# # selection_coef_long <- gather(selection_coef, scenario_population, value, factor_key=TRUE)
# # tapply(selection_coef_long$value, selection_coef_long$scenario_population, mean)
# # selection_coef_long_split = data.frame(str_split_fixed(selection_coef_long$scenario_population, fixed("."), 2))
# # selection_coef_long$scenario = selection_coef_long_split$X1 
# # selection_coef_long$population = selection_coef_long_split$X2 
# # selection_coef_long$coeff = "10 <"

# # strong_ben = selection_coef_long


# # dt = rbind(strong_del, weak_del, neu_del, neu_ben, weak_ben, strong_ben)
dt = rbind(strong_del, weak_del, neu, ben)

dt$scenario <- factor(dt$scenario, levels = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11", "X12", "X13", "X14", "X15", "X16", "X17", "X18", "X19", "X20", "X21", "X22", "X23", "X24"), labels = c("Mig=0%, Sb=0, pb=0%, Change=0%", "Mig=0%, Sb=0, pb=0%, Change=5%", "Mig=0%, Sb=0, pb=0%, Change=25%", "Mig=1%, Sb=0, pb=0%, Change=0%", "Mig=1%, Sb=0, pb=0%, Change=5%", "Mig=1%, Sb=0, pb=0%, Change=25%", "Mig=0%, Sb=1, pb=10%, Change=0%", "Mig=0%, Sb=1, pb=10%, Change=5%", "Mig=0%, Sb=1, pb=10%, Change=25%", "Mig=1%, Sb=1, pb=10%, Change=0%", "Mig=1%, Sb=1, pb=10%, Change=5%", "Mig=1%, Sb=1, pb=10%, Change=25%", "Mig=0%, Sb=10, pb=1%, Change=0%", "Mig=0%, Sb=10, pb=1%, Change=5%", "Mig=0%, Sb=10, pb=1%, Change=25%", "Mig=1%, Sb=10, pb=1%, Change=0%", "Mig=1%, Sb=10, pb=1%, Change=5%", "Mig=1%, Sb=10, pb=1%, Change=25%", "Mig=0%, Sb=100, pb=0.1%, Change=0%", "Mig=0%, Sb=100, pb=0.1%, Change=5%", "Mig=0%, Sb=100, pb=0.1%, Change=25%", "Mig=1%, Sb=100, pb=0.1%, Change=0%", "Mig=1%, Sb=100, pb=0.1%, Change=5%", "Mig=1%, Sb=100, pb=0.1%, Change=25%"), ordered = TRUE)

mean_boot = data.frame(tapply(dt$value, list(dt$scenario_population, dt$coeff), mean, na.rm=T))
CI5_boot  = data.frame(tapply(dt$value, list(dt$scenario_population, dt$coeff), quantile, p=c(0.05), na.rm=T))
CI95_boot  = data.frame(tapply(dt$value, list(dt$scenario_population, dt$coeff), quantile, p=c(0.95), na.rm=T))

mean_boot$s_p = rownames(mean_boot)
mean_boot = mean_boot %>%
  separate(s_p, c("scenario", "population"))
mean_boot_long = gather(mean_boot, coeff, value, -c(scenario, population), factor_key=TRUE)

CI5_boot$s_p = rownames(CI5_boot)
CI5_boot = CI5_boot %>%
  separate(s_p, c("scenario", "population"))
CI5_boot_long = gather(CI5_boot, coeff, CI5, -c(scenario, population), factor_key=TRUE)

CI95_boot$s_p = rownames(CI95_boot)
CI95_boot = CI95_boot %>%
  separate(s_p, c("scenario", "population"))
CI95_boot_long = gather(CI95_boot, coeff, CI95, -c(scenario, population), factor_key=TRUE)

dt_plot = inner_join(mean_boot_long, CI5_boot_long)
dt_plot = inner_join(dt_plot, CI95_boot_long)


dt_plot$coeff = gsub("ben", "> 1", dt_plot$coeff)
dt_plot$coeff = gsub("neu", "(-1, 1]", dt_plot$coeff)
dt_plot$coeff = gsub("weakly_del", "(-10, -1]", dt_plot$coeff)
dt_plot$coeff = gsub("strong_del", "< -10", dt_plot$coeff)
dt_plot$coeff <- factor(dt_plot$coeff, levels = c("< -10", "(-10, -1]", "(-1, 1]", "> 1"))

#dt_plot$scenario <- factor(dt_plot$scenario, levels = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11", "X12", "X13", "X14", "X15", "X16", "X17", "X18"), labels = c("Mig=0%, Sb=1, pb=10%, Change=0%", "Mig=0%, Sb=1, pb=10%, Change=5%", "Mig=0%, Sb=1, pb=10%, Change=25%", "Mig=1%, Sb=1, pb=10%, Change=0%", "Mig=1%, Sb=1, pb=10%, Change=5%", "Mig=1%, Sb=1, pb=10%, Change=25%", "Mig=0%, Sb=10, pb=1%, Change=0%", "Mig=0%, Sb=10, pb=1%, Change=5%", "Mig=0%, Sb=10, pb=1%, Change=25%", "Mig=1%, Sb=10, pb=1%, Change=0%", "Mig=1%, Sb=10, pb=1%, Change=5%", "Mig=1%, Sb=10, pb=1%, Change=25%", "Mig=0%, Sb=100, pb=0.1%, Change=0%", "Mig=0%, Sb=100, pb=0.1%, Change=5%", "Mig=0%, Sb=100, pb=0.1%, Change=25%", "Mig=1%, Sb=100, pb=0.1%, Change=0%", "Mig=1%, Sb=100, pb=0.1%, Change=5%", "Mig=1%, Sb=100, pb=0.1%, Change=25%"), ordered = TRUE)
dt_plot$scenario <- factor(dt_plot$scenario, levels = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11", "X12", "X13", "X14", "X15", "X16", "X17", "X18", "X19", "X20", "X21", "X22", "X23", "X24"), labels = c("Mig=0%, Sb=0, pb=0%, Change=0%", "Mig=0%, Sb=0, pb=0%, Change=5%", "Mig=0%, Sb=0, pb=0%, Change=25%", "Mig=1%, Sb=0, pb=0%, Change=0%", "Mig=1%, Sb=0, pb=0%, Change=5%", "Mig=1%, Sb=0, pb=0%, Change=25%", "Mig=0%, Sb=1, pb=10%, Change=0%", "Mig=0%, Sb=1, pb=10%, Change=5%", "Mig=0%, Sb=1, pb=10%, Change=25%", "Mig=1%, Sb=1, pb=10%, Change=0%", "Mig=1%, Sb=1, pb=10%, Change=5%", "Mig=1%, Sb=1, pb=10%, Change=25%", "Mig=0%, Sb=10, pb=1%, Change=0%", "Mig=0%, Sb=10, pb=1%, Change=5%", "Mig=0%, Sb=10, pb=1%, Change=25%", "Mig=1%, Sb=10, pb=1%, Change=0%", "Mig=1%, Sb=10, pb=1%, Change=5%", "Mig=1%, Sb=10, pb=1%, Change=25%", "Mig=0%, Sb=100, pb=0.1%, Change=0%", "Mig=0%, Sb=100, pb=0.1%, Change=5%", "Mig=0%, Sb=100, pb=0.1%, Change=25%", "Mig=1%, Sb=100, pb=0.1%, Change=0%", "Mig=1%, Sb=100, pb=0.1%, Change=5%", "Mig=1%, Sb=100, pb=0.1%, Change=25%"), ordered = TRUE)

scaling = read.delim(file="../results/data_frames/polyDFE_diversity_Ne_scaling.txt", header = T)

scaling = scaling %>%
  separate(params, c("mig", "possel", "pchange"), "-")

expectation = data.frame()
for (POP in c("Domesticated", "Wild"))
{
  # POP = "Domesticated"
  for (i in seq(24))
  {
    # i = 2
    ROW = slice(scaling, i)
    POSSEL = ROW$possel
    
    if (POP == "Domesticated") 
    {
      SCALING = ROW$Domesticated/5000
    }
    
    if (POP == "Wild") 
    {
      SCALING = ROW$Wild/5000
    }
    
    if (POSSEL == 0)
    {
      Sb = 0
      pb = 0
    }
    
    if (POSSEL == 2)
    {
      Sb = 1
      pb = 0.1
    }
    
    if (POSSEL == 20)
    {
      Sb = 10
      pb = 0.01
    }
    
    if (POSSEL == 200)
    {
      Sb = 100
      pb = 0.001
    }
    
    neu = pgamma(q=c(1), shape=0.3, scale=(100/0.3)*SCALING)
    strong_del = 1-pgamma(q=c(10), shape=0.3, scale=(100/0.3)*SCALING)
    weak_del = 1 - neu - strong_del
    neu2 = pexp(q=c(1),rate=1/(Sb*SCALING))
    ben = 1-neu2
    
    BENEFICIAL = pb*ben
    NEUTRAL = pb*neu2 + (1-pb)*neu
    WEAKLY_DEL = (1-pb)*weak_del
    STRONG_DEL = (1-pb)*strong_del
    
    BENEFICIAL+NEUTRAL+WEAKLY_DEL+STRONG_DEL
    
    exp_1 = c(paste0(ROW$mig,"-",ROW$possel,"-",ROW$pchange), POP, "< -10", STRONG_DEL, STRONG_DEL, STRONG_DEL)
    exp_2 = c(paste0(ROW$mig,"-",ROW$possel,"-",ROW$pchange), POP, "(-10, -1]", WEAKLY_DEL, WEAKLY_DEL, WEAKLY_DEL)
    exp_3 = c(paste0(ROW$mig,"-",ROW$possel,"-",ROW$pchange), POP, "(-1, 1]", NEUTRAL, NEUTRAL, NEUTRAL)
    exp_4 = c(paste0(ROW$mig,"-",ROW$possel,"-",ROW$pchange), POP, "> 1", BENEFICIAL, BENEFICIAL, BENEFICIAL)
    
    results = rbind(exp_1, exp_2, exp_3, exp_4)
    expectation = rbind(expectation, results)
  }
}

expectation$V1 <- factor(expectation$V1, levels = c('0-0-0', '0-0-0.05', '0-0-0.25', '0.01-0-0', '0.01-0-0.05', '0.01-0-0.25', '0-2-0', '0-2-0.05', '0-2-0.25', '0.01-2-0', '0.01-2-0.05', '0.01-2-0.25', '0-20-0', '0-20-0.05', '0-20-0.25', '0.01-20-0', '0.01-20-0.05', '0.01-20-0.25', '0-200-0', '0-200-0.05', '0-200-0.25', '0.01-200-0', '0.01-200-0.05', '0.01-200-0.25'), labels = c("Mig=0%, Sb=0, pb=0%, Change=0%", "Mig=0%, Sb=0, pb=0%, Change=5%", "Mig=0%, Sb=0, pb=0%, Change=25%", "Mig=1%, Sb=0, pb=0%, Change=0%", "Mig=1%, Sb=0, pb=0%, Change=5%", "Mig=1%, Sb=0, pb=0%, Change=25%", "Mig=0%, Sb=1, pb=10%, Change=0%", "Mig=0%, Sb=1, pb=10%, Change=5%", "Mig=0%, Sb=1, pb=10%, Change=25%", "Mig=1%, Sb=1, pb=10%, Change=0%", "Mig=1%, Sb=1, pb=10%, Change=5%", "Mig=1%, Sb=1, pb=10%, Change=25%", "Mig=0%, Sb=10, pb=1%, Change=0%", "Mig=0%, Sb=10, pb=1%, Change=5%", "Mig=0%, Sb=10, pb=1%, Change=25%", "Mig=1%, Sb=10, pb=1%, Change=0%", "Mig=1%, Sb=10, pb=1%, Change=5%", "Mig=1%, Sb=10, pb=1%, Change=25%", "Mig=0%, Sb=100, pb=0.1%, Change=0%", "Mig=0%, Sb=100, pb=0.1%, Change=5%", "Mig=0%, Sb=100, pb=0.1%, Change=25%", "Mig=1%, Sb=100, pb=0.1%, Change=0%", "Mig=1%, Sb=100, pb=0.1%, Change=5%", "Mig=1%, Sb=100, pb=0.1%, Change=25%"), ordered = TRUE)

expectation$V2 = gsub("Domesticated", "Domestic", expectation$V2)
colnames(expectation) <- c("scenario", "population", "coeff", "value", "CI5", "CI95")
expectation$Value = "True"
dt_plot$Value = "Inferred"

dt_plot = rbind(dt_plot, expectation)
dt_plot$value = as.numeric(dt_plot$value)
dt_plot$CI5 = as.numeric(dt_plot$CI5)
dt_plot$CI95 = as.numeric(dt_plot$CI95)
dt_plot$Population = dt_plot$population

options(scipen=10000)
dfe_plot = ggplot(dt_plot, aes(y=as.numeric(value), x=as.factor(coeff), shape=Population, color=Value)) +
  geom_point(position=position_dodge(width = .5),size=4) +
  # geom_errorbar(aes(ymin = CI5, ymax = CI95)) + 
  geom_pointrange(aes(ymin=CI5, ymax=CI95), position=position_dodge(width = .5)) +
  scale_color_grey() +
  xlab("4Nes") +
  ylab("Fraction") +
  # scale_y_log10(breaks=c(0.0001, 0.001, 0.01, 0.1, 1)) +
  theme_bw() +  
  facet_wrap(~scenario, ncol = 3) +
  theme(text = element_text(size=20), axis.title.y=element_blank(), legend.position = "bottom")

ggsave(plot = dfe_plot, "../plots/polyDFE/AIC_weighted_DFE.png", width=20, height=15, dpi="retina")

tapply(dt_plot$value, list(paste0(dt_plot$scenario,"-",dt_plot$Value), dt_plot$population), sum, na.rm=T) 

