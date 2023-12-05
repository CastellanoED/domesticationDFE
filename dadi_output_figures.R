library("data.table")
library("dplyr")
library("tidyr")
library("ggplot2")
library(gridExtra)
library(openxlsx)
#library(ggpubr)

setwd("/xdisk/rgutenk/dcastellano/Projects/simulations/test6/scripts")

### Results using the original fs - no replicates, no CI ###
dt = fread("../results/data_frames/vourlaki_DFE_original_fs_best_fits.txt", fill = T)
colnames(dt) <- c("params", "loglk", "shape", "scale", "pb_wild", "Sb", "pchange", "pb_domesticated", "missid", "theta_nsyn")

p <- ggplot(dt, aes(y=params, x=loglk)) + 
  geom_boxplot() +
  theme_linedraw() #+

p

dt = dt %>%
  separate(params, c("mig", "true_possel", "true_pchange", "assumed_possel", "assumed_pb_num"), "-")

best_loglk = data.frame(tapply(dt$loglk, paste0(dt$mig,"-",dt$true_possel,"-",dt$true_pchange), max))
best_loglk$params = row.names(best_loglk)
colnames(best_loglk) <- c("loglk", "params")

dt$params = paste0(dt$mig,"-",dt$true_possel,"-",dt$true_pchange)
best_loglk$loglk = as.numeric(as.character(best_loglk$loglk))

dt_best_loglk = inner_join(best_loglk, dt)

p <- ggplot(dt_best_loglk, aes(y=params, x=pchange)) + 
  geom_boxplot() +
  theme_linedraw() #+
p

p <- ggplot(dt_best_loglk, aes(y=params, x=shape*scale)) + 
  geom_boxplot() +
  theme_linedraw() #+

p

p <- ggplot(dt_best_loglk, aes(y=params, x=shape)) + 
  geom_boxplot() +
  theme_linedraw() #+

p

p <- ggplot(dt_best_loglk, aes(y=params, x=pb_wild)) + 
  geom_boxplot() +
  theme_linedraw() #+

p


### Results using the bootstraped fs - with replicates and CI ###
dt = fread("../results/data_frames/vourlaki_DFE_bootstra_fs_best_fits.txt", fill = T)
dt = na.omit(dt) # removing the model with a single pb for both populations

colnames(dt) <- c("params", "loglk", "shape", "scale", "pb_wild", "Sb", "pchange", "pb_domesticated", "missid", "theta_nsyn")

dt = dt %>%
  separate(params, c("replicate", "mig", "true_possel", "true_pchange", "assumed_possel", "assumed_pb_num"), "-")

best_loglk = data.frame(tapply(dt$loglk, paste0(dt$replicate,"-",dt$mig,"-",dt$true_possel,"-",dt$true_pchange), max))
best_loglk$params = row.names(best_loglk)
colnames(best_loglk) <- c("loglk", "params")
best_loglk = best_loglk %>%
  separate(params, c("replicate", "mig", "true_possel", "true_pchange"), "-")

best_loglk$params = paste0(best_loglk$mig,"-",best_loglk$true_possel,"-",best_loglk$true_pchange)
dt$params = paste0(dt$mig,"-",dt$true_possel,"-",dt$true_pchange)

dt_best_loglk = inner_join(best_loglk, dt)
dt_best_loglk$params <- factor(dt_best_loglk$params, levels = c('0-2-0','0-2-0.05', '0-2-0.25', '0.01-2-0','0.01-2-0.05', '0.01-2-0.25', '0-20-0','0-20-0.05', '0-20-0.25', '0.01-20-0','0.01-20-0.05', '0.01-20-0.25', '0-200-0','0-200-0.05', '0-200-0.25', '0.01-200-0','0.01-200-0.05', '0.01-200-0.25'), labels = c("ID = 1", "ID = 2", "ID = 3", "ID = 4", "ID = 5", "ID = 6", "ID = 7", "ID = 8", "ID = 9", "ID = 10", "ID = 11", "ID = 12", "ID = 13", "ID = 14", "ID = 15", "ID = 16", "ID = 17", "ID = 18"), ordered = TRUE)
dt_best_loglk$ancestral_Ne <- dt_best_loglk$theta_nsyn / 80

tapply(dt_best_loglk$ancestral_Ne, dt_best_loglk$params, mean)
table(list(dt_best_loglk$true_possel, dt_best_loglk$assumed_possel, dt_best_loglk$true_pchange, dt_best_loglk$mig))
tapply(dt_best_loglk$pchange, dt_best_loglk$params, quantile, probs=c(0.025,0.975))

# mu*L is 20 for nsyn
(100 * 2/3  * 120   * 10000    * 2.5e-7)  * 4 * 5000 #expected theta nsyn under neutrality
#(it * nsyn * #loci * bp/locus * mu/site)  * heritable units * "2"  * N
#The "2" comes from the fact that two sequences that have diverged for time t are different by 2 * mu * t mutations, since both diverging lineages accumulate mutations.

### Comparison of Ne inference across populations and  methods ###
demog = fread("../results/data_frames/flexible_dem_best_fits.txt", fill = T)
colnames(demog) <- c("params", "loglk", "Tpre", "nuPre", "Tdiv", "nu1div", "nu2div", "T1F", "T2F", "nu1F", "nu2F", "mw2d", "md2w", "missid","theta_syn", "ratio_nsyn_syn")
demog$Wild = (demog$nu1F*demog$theta_syn) / 40
demog$Domesticated = (demog$nu2F*demog$theta_syn) / 40
demog$Ancestral = (demog$theta_syn) / 40
demog_Ne = dplyr::select(demog, "params", "Ancestral", "Wild", "Domesticated")
demog_Ne$params <- factor(demog_Ne$params, levels = c('0-2-0','0-2-0.05', '0-2-0.25', '0.01-2-0','0.01-2-0.05', '0.01-2-0.25', '0-20-0','0-20-0.05', '0-20-0.25', '0.01-20-0','0.01-20-0.05', '0.01-20-0.25', '0-200-0','0-200-0.05', '0-200-0.25', '0.01-200-0','0.01-200-0.05', '0.01-200-0.25'), labels = c("Mig=0%, Sb=1, pb=10%, Change=0%", "Mig=0%, Sb=1, pb=10%, Change=5%", "Mig=0%, Sb=1, pb=10%, Change=25%", "Mig=1%, Sb=1, pb=10%, Change=0%", "Mig=1%, Sb=1, pb=10%, Change=5%", "Mig=1%, Sb=1, pb=10%, Change=25%", "Mig=0%, Sb=10, pb=1%, Change=0%", "Mig=0%, Sb=10, pb=1%, Change=5%", "Mig=0%, Sb=10, pb=1%, Change=25%", "Mig=1%, Sb=10, pb=1%, Change=0%", "Mig=1%, Sb=10, pb=1%, Change=5%", "Mig=1%, Sb=10, pb=1%, Change=25%", "Mig=0%, Sb=100, pb=0.1%, Change=0%", "Mig=0%, Sb=100, pb=0.1%, Change=5%", "Mig=0%, Sb=100, pb=0.1%, Change=25%", "Mig=1%, Sb=100, pb=0.1%, Change=0%", "Mig=1%, Sb=100, pb=0.1%, Change=5%", "Mig=1%, Sb=100, pb=0.1%, Change=25%"), ordered = TRUE)

tapply(demog_Ne$Ancestral, demog_Ne$params, mean)
demog_Ne <- gather(demog_Ne, Population, Ne, Ancestral:Domesticated, factor_key=TRUE)
demog_Ne$Method = "dadi"
demog_Ne_polyDFE = fread("../results/data_frames/polyDFE_Ne.txt", header = T)
demog_Ne = rbind(demog_Ne, demog_Ne_polyDFE)

myplot = ggplot(demog_Ne, aes(x=Ne/5000, y = 1, color=Population, shape=Method)) +
  ggtitle("Ne dadi's best demographic model vs polyDFE") +
  geom_point(size=5,alpha=0.75) +
  scale_color_brewer(palette="Dark2") +
  facet_wrap(~ params, nrow = 6) +
  #scale_fill_grey() +
  geom_vline(xintercept = 1, linetype="dashed") +
  xlim(0,2) +
  theme_classic() +
  xlab("Relative Ne") +
  theme(text = element_text(size=15), axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position="bottom")
ggsave(myplot, file="../plots/dadi/Na.png", width=15.7, height=12, dpi=100)


### Plot dadi's DFE parameters ###

### pchange
A = ggplot(dt_best_loglk, aes(x=pchange)) +
  ggtitle("pc bootstrap distributions") +
  geom_histogram(position = "identity", bins=30) +
  geom_density(adjust=0.5, color="red") +
  #geom_boxplot() + 
  facet_wrap(~ params, nrow = 6) +
  #scale_fill_grey() +
  geom_vline(xintercept = 0.25, linetype="dashed") +
  geom_vline(xintercept = 0.05, linetype="dotdash") +
  geom_vline(xintercept = 0, linetype="dotted") +
  xlim(-0.1,0.5) +
  theme_classic() +
  xlab("pc") +
  theme(text = element_text(size=15), axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position="bottom")
A
ggsave(A, file="../plots/dadi/pchange.png", width=12, height=15.7, dpi=100)


dt_best_loglk <- dt_best_loglk %>% 
  mutate(scenario = case_when(true_possel == 2 ~   "Pervasive and nearly neutral",
                              true_possel == 20 ~  "Common and weak",
                              true_possel == 200 ~ "Rare and strong"))
dt_best_loglk$scenario <- factor(dt_best_loglk$scenario, levels = c("Pervasive and nearly neutral", "Common and weak", "Rare and strong"), ordered = TRUE)

dummy = fread("../plots/dummy_fig4.txt")
dummy$true_pchange = as.factor(dummy$true_pchange)
dummy$mig = as.factor(dummy$mig)
dummy$value = "True"

dt_fig4 = dplyr::select(dt_best_loglk, scenario, true_pchange, mig, pchange)
dt_fig4$value = "Inferred"
dt_fig4 = rbind(dt_fig4, dummy)

B = ggplot() +
  geom_violin(data=subset(dt_fig4, value == "Inferred"), aes(x = true_pchange, y = as.numeric(pchange), fill =mig)) +
  geom_boxplot(data=subset(dt_fig4, value == "True"), aes(x = true_pchange, y = as.numeric(pchange), color="black")) +
  theme_bw() +
  scale_fill_brewer(palette = "Greens") +
  #scale_color_brewer(palette = "Greens") +
  
  facet_wrap(~scenario, nrow=3) +
  #geom_hline(yintercept = 0.25, linetype="dashed") +
  #geom_hline(yintercept = 0.05, linetype="dotdash") +
  ylab("Inferred pc") + xlab("Simulated pc") + theme(text = element_text(size=25), legend.position = "none")
B
ggsave(B, file="../plots/dadi/Figure4.png", width=12, height=15.7, dpi=100)

ggplot(dt_best_loglk, aes(x=pchange)) +
  #ggtitle("pc bootstrap distributions") +
  #geom_point(adjust=2, alpha=0.5) +
  geom_histogram(position = "identity", alpha = 0.5, bins=30) +
  #geom_density(alpha = 0.5) + 
  
  facet_wrap(~ params, nrow = 6) +
  #scale_fill_brewer(palette = "Blues") +
  #scale_fill_grey() +
  geom_vline(xintercept = 0.25, linetype="dashed") +
  geom_vline(xintercept = 0.05, linetype="dotdash") +
  geom_vline(xintercept = 0, linetype="dotted") +
  xlim(0,0.5) +
  theme_classic() +
  ylab("true pc") + xlab("inferred pc")

ggplot(dt_best_loglk, aes(y=params, x=pchange, fill=true_pchange ) ) + 
  geom_boxplot() +
  scale_fill_brewer(palette="BuPu") + xlim(0,.5) + 
  geom_vline(xintercept = 0.25, linetype="dashed") +
  geom_vline(xintercept = 0.05, linetype="dotdash") +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme_linedraw() #+

### shape
myplot = ggplot(dt_best_loglk, aes(x=shape)) +
  ggtitle("B") +
  geom_density(adjust=2, alpha=0.5) +
  facet_wrap(~ params, nrow = 6) +
  #scale_fill_grey() +
  geom_vline(xintercept = 0.3, linetype="dashed") +
  xlim(0,1.6) +
  theme_classic() +
  xlab("") +
  theme(text = element_text(size=25), axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position="bottom")
panelB = myplot
ggsave(myplot, file="../plots/dadi/b.png", width=12, height=15.7, dpi=100)

### Sb
myplot = ggplot(dt_best_loglk, aes(x=log10(Sb))) +
  ggtitle("Sb bootstrap distributions") +
  geom_density(adjust=2, alpha=0.5) +
  facet_wrap(~ params, nrow = 6) +
  #scale_fill_grey() +
  geom_vline(xintercept = log10(1), linetype="dashed") +
  geom_vline(xintercept = log10(10), linetype="dashed") +
  geom_vline(xintercept = log10(100), linetype="dashed") +
  #xlim(0,1.6) +
  theme_classic() +
  xlab("log10(Sb)") +
  theme(text = element_text(size=15), axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position="bottom")
ggsave(myplot, file="../plots/dadi/Sb.png", width=12, height=15.7, dpi=100)

### sd
myplot = ggplot(dt_best_loglk, aes(x=2*20*((shape*scale)/theta_nsyn)   )) +
  ggtitle("D") +
  geom_density(adjust=2, alpha=0.5) +
  facet_wrap(~ params, nrow = 6) +
  #scale_fill_grey() +
  geom_vline(xintercept = 0.01, linetype="dashed") +
  xlim(0,0.02) +
  theme_classic() +
  xlab("") +
  theme(text = element_text(size=25), axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position="bottom", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
panelD = myplot
ggsave(myplot, file="../plots/dadi/sd.png", width=12, height=15.7, dpi=100)

### pb
dt_pbs_wide = dplyr::select(dt_best_loglk, params, pb_wild, pb_domesticated, pchange)
hist(dt_pbs_wide$pb_wild)
dt_pbs_wide$pb_domesticated = dt_pbs_wide$pb_wild*(1-dt_pbs_wide$pchange)+dt_pbs_wide$pb_domesticated*dt_pbs_wide$pchange
hist(dt_pbs_wide$pb_wild)
dt_pbs_long <- gather(dt_pbs_wide, Population, pb, pb_wild:pb_domesticated, factor_key=TRUE)
dt_pbs_long$Population = gsub("pb_wild", "Wild", dt_pbs_long$Population)
dt_pbs_long$Population = gsub("pb_domesticated", "Domesticated", dt_pbs_long$Population)

myplot = ggplot(dt_pbs_long, aes(x=log10(pb+0.001), fill=Population)) +
  ggtitle("pb bootstrap distributions") +
  #geom_density(adjust=2, alpha=0.5) +
  geom_boxplot( alpha=0.5) +
  facet_wrap(~ params, nrow = 6) +
  scale_fill_grey() +
  geom_vline(xintercept = -1, linetype="dashed") +
  geom_vline(xintercept = -2, linetype="dashed") +
  geom_vline(xintercept = -3, linetype="dashed") +
  #xlim(0,1.6) +
  theme_classic() +
  xlab("log10(pb)") +
  theme(text = element_text(size=15), axis.text.y=element_blank(), axis.title.y=element_blank(),legend.position="bottom")
ggsave(myplot, file="../plots/dadi/pb.png", width=12, height=15.7, dpi=100)


### Calculating Pn/Ps across bootstrap replicates

### Estimates for each independent inference unit ###
pNpS_all_sims = data.frame()

for (it in seq(0, 99))
{
  for (mig in c(0, 0.01))
  {
    for (possel in c(2, 20, 200))
    {
      for (pchange in c(0, 0.05, 0.25))
      {
        dt = fread(paste0("../results/SFS/by_dadi/Bootstraps/sim_", mig,"-", possel, "-", pchange, "-nsyn_boot.", it, ".domesticated_subsampled.1d.fs"), skip = 1, nrows = 1)
        dt$sum = rowSums(dt)
        dt$sum = dt$sum  - (dt$V1 + dt$V21)
        pN = dt$sum
        
        dt = fread(paste0("../results/SFS/by_dadi/Bootstraps/sim_", mig,"-", possel, "-", pchange, "-syn_boot.", it, ".domesticated_subsampled.1d.fs"), skip = 1, nrows = 1)
        dt$sum = rowSums(dt)
        dt$sum = dt$sum  - (dt$V1 + dt$V21)
        pS = dt$sum
        pop = "domesticated"
        
        pNpS = data.frame(it, mig, possel, pchange, pop, pN, pS)
        pNpS_all_sims = rbind(pNpS_all_sims, pNpS)
        
        dt = fread(paste0("../results/SFS/by_dadi/Bootstraps/sim_", mig,"-", possel, "-", pchange,  "-nsyn_boot.", it, ".wild_subsampled.1d.fs"), skip = 1, nrows = 1)
        dt$sum = rowSums(dt)
        dt$sum = dt$sum  - (dt$V1 + dt$V21)
        pN = dt$sum
        
        dt = fread(paste0("../results/SFS/by_dadi/Bootstraps/sim_", mig,"-", possel, "-", pchange, "-syn_boot.", it, ".wild_subsampled.1d.fs"), skip = 1, nrows = 1)
        dt$sum = rowSums(dt)
        dt$sum = dt$sum  - (dt$V1 + dt$V21)
        pS = dt$sum
        pop = "wild"
        
        pNpS = data.frame(it, mig, possel, pchange, pop, pN, pS)
        pNpS_all_sims = rbind(pNpS_all_sims, pNpS)
        
      }
    }
  }
}

pNpS_all_sims$ratio = pNpS_all_sims$pN/(pNpS_all_sims$pS*2)
hist(subset(pNpS_all_sims, pop=="wild")$ratio)
hist(subset(pNpS_all_sims, pop=="domesticated")$ratio)

hist(subset(pNpS_all_sims, pop=="wild")$pS)
hist(subset(pNpS_all_sims, pop=="domesticated")$pS)

pNpS_all_sims$params = paste0(pNpS_all_sims$mig, "-", pNpS_all_sims$possel, "-", pNpS_all_sims$pchange)
pNpS_all_sims$params <- factor(pNpS_all_sims$params, levels = c('0-2-0','0.01-2-0', '0-20-0','0.01-20-0', '0-200-0','0.01-200-0', '0-2-0.05','0.01-2-0.05', '0-20-0.05','0.01-20-0.05', '0-200-0.05','0.01-200-0.05', '0-2-0.25','0.01-2-0.25', '0-20-0.25','0.01-20-0.25', '0-200-0.25','0.01-200-0.25'),ordered = TRUE)

p <- ggplot(pNpS_all_sims, aes(y=params, x=ratio, fill=pop ) ) + 
  geom_boxplot() +
  scale_fill_brewer(palette="Dark2") +
  theme_linedraw() + xlab("Pn/Ps")
p

p <- ggplot(pNpS_all_sims, aes(y=params, x=pS, fill=pop ) ) + 
  geom_boxplot() +
  scale_fill_brewer(palette="Dark2") +
  theme_linedraw() + xlab("Ps")
p


### Estimates for each independent simulation run ###
pNpS_all_sims = data.frame()

for (it in seq(1, 100))
{
  for (mig in c(0, 0.01))
  {
    for (possel in c(2, 20, 200))
    {
      for (pchange in c(0, 0.05, 0.25))
      {
        dt_nsyn = fread(paste0("../results/VCF/sim_", mig,"-", possel, "-", pchange, "-", it, "_nsyn.vcf"))
        dt_syn =  fread(paste0("../results/VCF/sim_", mig,"-", possel, "-", pchange, "-", it, "_syn.vcf"))
        
        dt_dom = dplyr::select(dt_nsyn, POS, D1:D9)
        dt_dom$sum = rowSums(dt_dom)
        dt_dom$sum = dt_dom$sum  - dt_dom$POS
        dt_dom = subset(dt_dom, sum > 0 & sum < 40)
        pN = nrow(dt_dom)

        dt_dom = dplyr::select(dt_syn, POS, D1:D9)
        dt_dom$sum = rowSums(dt_dom)
        dt_dom$sum = dt_dom$sum  - dt_dom$POS
        dt_dom = subset(dt_dom, sum > 0 & sum < 40)
        pS = nrow(dt_dom)
        pop = "Domestic"
        
        pNpS = data.frame(it, mig, possel, pchange, pop, pN, pS)
        pNpS_all_sims = rbind(pNpS_all_sims, pNpS)
        
        dt_wild = dplyr::select(dt_nsyn, POS, W1:W9)
        dt_wild$sum = rowSums(dt_wild)
        dt_wild$sum = dt_wild$sum  - dt_wild$POS
        dt_wild = subset(dt_wild, sum > 0 & sum < 40)
        pN = nrow(dt_wild)
        
        dt_wild = dplyr::select(dt_syn, POS, W1:W9)
        dt_wild$sum = rowSums(dt_wild)
        dt_wild$sum = dt_wild$sum  - dt_wild$POS
        dt_wild = subset(dt_wild, sum > 0 & sum < 40)
        pS = nrow(dt_wild)
        pop = "Wild"
        
        pNpS = data.frame(it, mig, possel, pchange, pop, pN, pS)
        pNpS_all_sims = rbind(pNpS_all_sims, pNpS)
        
      }
    }
  }
}

pNpS_all_sims$ratio = pNpS_all_sims$pN/(pNpS_all_sims$pS*2)

pNpS_all_sims$params = paste0(pNpS_all_sims$mig, "-", pNpS_all_sims$possel, "-", pNpS_all_sims$pchange)
pNpS_all_sims$params <- factor(pNpS_all_sims$params, levels = c('0-2-0','0-2-0.05', '0-2-0.25', '0.01-2-0','0.01-2-0.05', '0.01-2-0.25', '0-20-0','0-20-0.05', '0-20-0.25', '0.01-20-0','0.01-20-0.05', '0.01-20-0.25', '0-200-0','0-200-0.05', '0-200-0.25', '0.01-200-0','0.01-200-0.05', '0.01-200-0.25'), labels = c("Mig=0%, Sb=1, pb=10%, Change=0%", "Mig=0%, Sb=1, pb=10%, Change=5%", "Mig=0%, Sb=1, pb=10%, Change=25%", "Mig=1%, Sb=1, pb=10%, Change=0%", "Mig=1%, Sb=1, pb=10%, Change=5%", "Mig=1%, Sb=1, pb=10%, Change=25%", "Mig=0%, Sb=10, pb=1%, Change=0%", "Mig=0%, Sb=10, pb=1%, Change=5%", "Mig=0%, Sb=10, pb=1%, Change=25%", "Mig=1%, Sb=10, pb=1%, Change=0%", "Mig=1%, Sb=10, pb=1%, Change=5%", "Mig=1%, Sb=10, pb=1%, Change=25%", "Mig=0%, Sb=100, pb=0.1%, Change=0%", "Mig=0%, Sb=100, pb=0.1%, Change=5%", "Mig=0%, Sb=100, pb=0.1%, Change=25%", "Mig=1%, Sb=100, pb=0.1%, Change=0%", "Mig=1%, Sb=100, pb=0.1%, Change=5%", "Mig=1%, Sb=100, pb=0.1%, Change=25%"), ordered = TRUE)

A <- ggplot(pNpS_all_sims, aes(x=ratio, fill=pop ) ) + 
  geom_histogram(alpha=0.5, position="identity") +
  #geom_boxplot() + 
  facet_wrap(~ params, nrow = 6) +
  scale_fill_brewer(palette = "Dark2") +
  #geom_vline(xintercept = 0.25, linetype="dashed") +
  #geom_vline(xintercept = 0.05, linetype="dotdash") +
  #geom_vline(xintercept = 0, linetype="dotted") +
  #xlim(0,1) +
  theme_bw() +
  xlab("Pn/Ps") +
  theme(text = element_text(size=15), axis.title.y=element_blank(),legend.position="none")
A

armonic_number = sum(1/seq(39))
L = 120/3 * 10000

B <- ggplot(pNpS_all_sims, aes(x=pS/L/armonic_number, fill=pop ) ) + 
  geom_histogram(alpha=0.5, position="identity") +
  #geom_boxplot() + 
  facet_wrap(~ params, nrow = 6) +
  scale_fill_brewer(palette = "Dark2") +
  #geom_vline(xintercept = 0.25, linetype="dashed") +
  #geom_vline(xintercept = 0.05, linetype="dotdash") +
  #geom_vline(xintercept = 0, linetype="dotted") +
  #xlim(0,1) +
  theme_bw() +
  xlab(expression(theta[syn])) +
  theme(text = element_text(size=15), axis.title.y=element_blank(),legend.position="none")
B

myplot = grid.arrange(B, A, ncol = 2)
ggsave(myplot, file="../plots/diversity_PnPs.png", width=24, height=15, dpi=200)


### Supplementary Fig1 ###
pNpS_all_sims = data.frame()

for (it in seq(1, 100))
{
  for (mig in c(0, 0.01))
  {
    for (possel in c(2, 20, 200))
    {
      for (pchange in c(0, 0.05, 0.25))
      {
        # it = 1
        # mig = 0
        # possel = 20
        # pchange = 0.25
        
        dt_sites = fread(paste0("../results/fitness_position_matrix/sim_", mig,"-", possel, "-", pchange, "-", it, "_fitness_position_matrix.txt"))
        dt_sites = dplyr::select(dt_sites, -"s.p1")
        colnames(dt_sites) <- c("POS", "TYPE", "WILD_s", "DOMESTIC_s")
        
        dt_nsyn = fread(paste0("../results/VCF/sim_", mig,"-", possel, "-", pchange, "-", it, "_nsyn.vcf"))
        dt_syn =  fread(paste0("../results/VCF/sim_", mig,"-", possel, "-", pchange, "-", it, "_syn.vcf"))
        
        dt_nsyn = inner_join(dt_nsyn, dt_sites)
        dt_syn  = inner_join( dt_syn, dt_sites)
        
        dt_dom = dplyr::select(dt_nsyn, POS, TYPE, D1:D9)
        dt_dom$sum = rowSums(dt_dom)
        dt_dom$sum = dt_dom$sum  - dt_dom$POS  - dt_dom$TYPE
        dt_dom = subset(dt_dom, sum > 0)
        pN = nrow(dt_dom)
        nsyn_SFS_dom = table(dt_dom$sum)
        
        dt_dom = dplyr::select(dt_syn, POS, TYPE, D1:D9)
        dt_dom$sum = rowSums(dt_dom)
        dt_dom$sum = dt_dom$sum  - dt_dom$POS  - dt_dom$TYPE
        dt_dom = subset(dt_dom, sum > 0)
        pS = nrow(dt_dom)
        pop = "Domestic"
        syn_SFS_dom = table(dt_dom$sum)
        
        pNpS = data.frame(it, mig, possel, pchange, pop, pN, pS)
        pNpS_all_sims = rbind(pNpS_all_sims, pNpS)
        
        dt_wild = dplyr::select(dt_nsyn, POS, W1:W9)
        dt_wild$sum = rowSums(dt_wild)
        dt_wild$sum = dt_wild$sum  - dt_wild$POS  - dt_wild$TYPE
        dt_wild = subset(dt_wild, sum > 0)
        pN = nrow(dt_wild)
        
        dt_wild = dplyr::select(dt_syn, POS, W1:W9)
        dt_wild$sum = rowSums(dt_wild)
        dt_wild$sum = dt_wild$sum  - dt_wild$POS  - dt_wild$TYPE
        dt_wild = subset(dt_wild, sum > 0)
        pS = nrow(dt_wild)
        pop = "Wild"
        
        pNpS = data.frame(it, mig, possel, pchange, pop, pN, pS)
        pNpS_all_sims = rbind(pNpS_all_sims, pNpS)
        
      }
    }
  }
}


### Discretized DFE
pop_size = data.frame(tapply(dt_best_loglk$theta_nsyn, paste0(dt_best_loglk$mig, "-", dt_best_loglk$true_possel, "-", dt_best_loglk$true_pchange), mean))
pop_size$params = row.names(pop_size)
colnames(pop_size) <- c("theta_nsyn", "params")
pop_size = pop_size %>%
  separate(params, c("mig", "possel", "pchange"), "-")

pop_size$scaling = pop_size$theta_nsyn/4e+05 
pop_size$scaling = 1

MEAN_SD = 0.01

POP = "BOTH"
expectation = data.frame()
  for (i in seq(18))
  {
    # i = 18
    ROW = slice(pop_size, i)
    POSSEL = ROW$possel
    
    SCALING = ROW$scaling
  
    if (POSSEL == 2)
    {
      Sb = MEAN_SD/100
      pb = 0.1
    }
    
    if (POSSEL == 20)
    {
      Sb = MEAN_SD/10
      pb = 0.01
    }
    
    if (POSSEL == 200)
    {
      Sb = MEAN_SD/1
      pb = 0.001
    }
    
    neu = pgamma(q=c(0.001), shape=0.3, scale=(MEAN_SD/0.3)*SCALING)
    strong_del = 1-pgamma(q=c(0.01), shape=0.3, scale=(MEAN_SD/0.3)*SCALING)
    weak_del = 1 - neu - strong_del
    neu2 = pexp(q=c(0),rate=1/(Sb*SCALING))
    ben = 1-neu2
    
    BENEFICIAL = pb*ben
    NEUTRAL = (1-pb)*neu
    WEAKLY_DEL = (1-pb)*weak_del
    STRONG_DEL = (1-pb)*strong_del
    
    BENEFICIAL+NEUTRAL+WEAKLY_DEL+STRONG_DEL
    
    exp_1 = c(paste0(ROW$mig,"-",ROW$possel,"-",ROW$pchange), POP, "< -0.01", STRONG_DEL, STRONG_DEL, STRONG_DEL, SCALING)
    exp_2 = c(paste0(ROW$mig,"-",ROW$possel,"-",ROW$pchange), POP, "(-0.01, -0.001]", WEAKLY_DEL, WEAKLY_DEL, WEAKLY_DEL, SCALING)
    exp_3 = c(paste0(ROW$mig,"-",ROW$possel,"-",ROW$pchange), POP, "(-0.001, 0]", NEUTRAL, NEUTRAL, NEUTRAL, SCALING)
    exp_4 = c(paste0(ROW$mig,"-",ROW$possel,"-",ROW$pchange), POP, "> 0", BENEFICIAL, BENEFICIAL, BENEFICIAL, SCALING)
    
    results = rbind(exp_1, exp_2, exp_3, exp_4)
    expectation = rbind(expectation, results)
  }


expectation$V1 <- factor(expectation$V1, levels = c('0-2-0', '0-2-0.05', '0-2-0.25', '0.01-2-0', '0.01-2-0.05', '0.01-2-0.25', '0-20-0', '0-20-0.05', '0-20-0.25', '0.01-20-0', '0.01-20-0.05', '0.01-20-0.25', '0-200-0', '0-200-0.05', '0-200-0.25', '0.01-200-0', '0.01-200-0.05', '0.01-200-0.25'), labels = c("Mig=0%, Sb=1, pb=10%, Change=0%", "Mig=0%, Sb=1, pb=10%, Change=5%", "Mig=0%, Sb=1, pb=10%, Change=25%", "Mig=1%, Sb=1, pb=10%, Change=0%", "Mig=1%, Sb=1, pb=10%, Change=5%", "Mig=1%, Sb=1, pb=10%, Change=25%", "Mig=0%, Sb=10, pb=1%, Change=0%", "Mig=0%, Sb=10, pb=1%, Change=5%", "Mig=0%, Sb=10, pb=1%, Change=25%", "Mig=1%, Sb=10, pb=1%, Change=0%", "Mig=1%, Sb=10, pb=1%, Change=5%", "Mig=1%, Sb=10, pb=1%, Change=25%", "Mig=0%, Sb=100, pb=0.1%, Change=0%", "Mig=0%, Sb=100, pb=0.1%, Change=5%", "Mig=0%, Sb=100, pb=0.1%, Change=25%", "Mig=1%, Sb=100, pb=0.1%, Change=0%", "Mig=1%, Sb=100, pb=0.1%, Change=5%", "Mig=1%, Sb=100, pb=0.1%, Change=25%"), ordered = TRUE)
colnames(expectation) <- c("scenario", "population", "coeff", "value", "CI5", "CI95", "scaling")
expectation$Value = "True"
expectation$value = as.numeric(expectation$value)
expectation$CI5 = as.numeric(expectation$CI5)
expectation$CI95 = as.numeric(expectation$CI95)
expectation$Population = expectation$population
expectation$coeff <- factor(expectation$coeff, levels = c("< -0.01", "(-0.01, -0.001]", "(-0.001, 0]", "> 0"), ordered = TRUE)

options(scipen=10000)
dfe_plot = ggplot(expectation, aes(y=as.numeric(value), x=as.factor(coeff), shape=Population, color=Value)) +
  geom_point(position=position_dodge(width = .5),size=3) +
  # geom_errorbar(aes(ymin = CI5, ymax = CI95)) + 
  geom_pointrange(aes(ymin=CI5, ymax=CI95), position=position_dodge(width = .5)) +
  scale_color_grey() +
  xlab("s") +
  ylab("Fraction") +
  ylim(0,0.5) +
  #scale_y_log10(breaks=c(0.0001, 0.001, 0.01, 0.1, 0.5, 1)) +
  theme_bw() +  
  facet_wrap(~scenario, nrow = 6, ncol = 3) +
  theme(text = element_text(size=20), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom")
ggsave(plot = dfe_plot, "../plots/polyDFE/true_discretized_DFE.png", width=20, height=15, dpi="retina")

tapply(expectation$value, list(paste0(expectation$scenario,"-",expectation$Value), expectation$population), sum, na.rm=T) 


### Comparing true vs inferred discretized DFE
dt_int = dt_best_loglk
dt_int$params <- factor(dt_int$params, levels = c("Mig=0%, Sb=1, pb=10%, Change=0%", "Mig=0%, Sb=1, pb=10%, Change=5%", "Mig=0%, Sb=1, pb=10%, Change=25%", "Mig=1%, Sb=1, pb=10%, Change=0%", "Mig=1%, Sb=1, pb=10%, Change=5%", "Mig=1%, Sb=1, pb=10%, Change=25%", "Mig=0%, Sb=10, pb=1%, Change=0%", "Mig=0%, Sb=10, pb=1%, Change=5%", "Mig=0%, Sb=10, pb=1%, Change=25%", "Mig=1%, Sb=10, pb=1%, Change=0%", "Mig=1%, Sb=10, pb=1%, Change=5%", "Mig=1%, Sb=10, pb=1%, Change=25%", "Mig=0%, Sb=100, pb=0.1%, Change=0%", "Mig=0%, Sb=100, pb=0.1%, Change=5%", "Mig=0%, Sb=100, pb=0.1%, Change=25%", "Mig=1%, Sb=100, pb=0.1%, Change=0%", "Mig=1%, Sb=100, pb=0.1%, Change=5%", "Mig=1%, Sb=100, pb=0.1%, Change=25%"), labels = c('0-2-0','0-2-0.05', '0-2-0.25', '0.01-2-0','0.01-2-0.05', '0.01-2-0.25', '0-20-0','0-20-0.05', '0-20-0.25', '0.01-20-0','0.01-20-0.05', '0.01-20-0.25', '0-200-0','0-200-0.05', '0-200-0.25', '0.01-200-0','0.01-200-0.05', '0.01-200-0.25'), ordered = TRUE)

dt_int$sd = 2*20*((dt_int$shape*dt_int$scale)/dt_int$theta_nsyn)
dt_int$neu = pgamma(q=c(0.001), shape=dt_int$shape, scale=(dt_int$sd/dt_int$shape))
dt_int$strong_del = 1-pgamma(q=c(0.01), shape=dt_int$shape, scale=(dt_int$sd/dt_int$shape))
dt_int$weak_del = 1 - dt_int$neu - dt_int$strong_del
dt_int$ben_dom = dt_int$pb_domesticated*dt_int$pchange+(1-dt_int$pchange)*dt_int$pb_wild
dt_int$ben_wild = dt_int$pb_wild

dt_int$BENEFICIAL_Domesticated = dt_int$ben_dom
dt_int$NEUTRAL_Domesticated = (1-dt_int$ben_dom)*dt_int$neu
dt_int$WEAKLY_Domesticated = (1-dt_int$ben_dom)*dt_int$weak_del
dt_int$STRONG_Domesticated = (1-dt_int$ben_dom)*dt_int$strong_del

dt_int$BENEFICIAL_Wild = dt_int$ben_wild
dt_int$NEUTRAL_Wild = (1-dt_int$ben_wild)*dt_int$neu
dt_int$WEAKLY_Wild = (1-dt_int$ben_wild)*dt_int$weak_del
dt_int$STRONG_Wild = (1-dt_int$ben_wild)*dt_int$strong_del

dt_int = dplyr::select(dt_int, params, BENEFICIAL_Domesticated:STRONG_Wild)
dt_int_long <- gather(dt_int, Coeff_Population, Value, BENEFICIAL_Domesticated:STRONG_Wild, factor_key=TRUE)
dt_int_long = dt_int_long %>%
  separate(Coeff_Population, c("coeff", "Population"), "_")

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
dt_int$scenario <- factor(dt_int$scenario, levels = c('0-2-0','0-2-0.05', '0-2-0.25', '0.01-2-0','0.01-2-0.05', '0.01-2-0.25', '0-20-0','0-20-0.05', '0-20-0.25', '0.01-20-0','0.01-20-0.05', '0.01-20-0.25', '0-200-0','0-200-0.05', '0-200-0.25', '0.01-200-0','0.01-200-0.05', '0.01-200-0.25'), labels = c("Mig=0%, Sb=1, pb=10%, Change=0%", "Mig=0%, Sb=1, pb=10%, Change=5%", "Mig=0%, Sb=1, pb=10%, Change=25%", "Mig=1%, Sb=1, pb=10%, Change=0%", "Mig=1%, Sb=1, pb=10%, Change=5%", "Mig=1%, Sb=1, pb=10%, Change=25%", "Mig=0%, Sb=10, pb=1%, Change=0%", "Mig=0%, Sb=10, pb=1%, Change=5%", "Mig=0%, Sb=10, pb=1%, Change=25%", "Mig=1%, Sb=10, pb=1%, Change=0%", "Mig=1%, Sb=10, pb=1%, Change=5%", "Mig=1%, Sb=10, pb=1%, Change=25%", "Mig=0%, Sb=100, pb=0.1%, Change=0%", "Mig=0%, Sb=100, pb=0.1%, Change=5%", "Mig=0%, Sb=100, pb=0.1%, Change=25%", "Mig=1%, Sb=100, pb=0.1%, Change=0%", "Mig=1%, Sb=100, pb=0.1%, Change=5%", "Mig=1%, Sb=100, pb=0.1%, Change=25%"), ordered = TRUE)

dt_int$scaling = 1
dt_int$Value = "Inferred"
dt_int$population = dt_int$Population

expectation = rbind(expectation, dt_int)
expectation$value = as.numeric(expectation$value)
expectation$CI5 = as.numeric(expectation$CI5)
expectation$CI95 = as.numeric(expectation$CI95)
expectation$coeff <- factor(expectation$coeff, levels = c("< -0.01", "(-0.01, -0.001]", "(-0.001, 0]", "> 0"), ordered = TRUE)

write.table(expectation, file="../results/data_frames/discretized_DFE_dadi_vs_true.txt", col.names = T, row.names = F, sep = "\t", quote = F, append = F)


dt_int_dadi = fread(file = "../results/data_frames/discretized_DFE_dadi_vs_true.txt", header = T)
dt_int_polyDFE = fread(file = "../results/data_frames/discretized_DFE_polyDFE.txt", header = T)
dt_int_dadi$Method = "dadi"

expectation = rbind(dt_int_dadi, dt_int_polyDFE)

expectation$Population = gsub("BOTH", "Both", expectation$Population)
expectation$Method[expectation$population == "BOTH" & expectation$Method == 'dadi'] <- "True"
expectation$coeff <- factor(expectation$coeff, levels = c("< -0.01", "(-0.01, -0.001]", "(-0.001, 0]", "> 0"), ordered = TRUE)
expectation$scenario <- factor(expectation$scenario, levels = c("Mig=0%, Sb=1, pb=10%, Change=0%", "Mig=0%, Sb=1, pb=10%, Change=5%", "Mig=0%, Sb=1, pb=10%, Change=25%", "Mig=1%, Sb=1, pb=10%, Change=0%", "Mig=1%, Sb=1, pb=10%, Change=5%", "Mig=1%, Sb=1, pb=10%, Change=25%", "Mig=0%, Sb=10, pb=1%, Change=0%", "Mig=0%, Sb=10, pb=1%, Change=5%", "Mig=0%, Sb=10, pb=1%, Change=25%", "Mig=1%, Sb=10, pb=1%, Change=0%", "Mig=1%, Sb=10, pb=1%, Change=5%", "Mig=1%, Sb=10, pb=1%, Change=25%", "Mig=0%, Sb=100, pb=0.1%, Change=0%", "Mig=0%, Sb=100, pb=0.1%, Change=5%", "Mig=0%, Sb=100, pb=0.1%, Change=25%", "Mig=1%, Sb=100, pb=0.1%, Change=0%", "Mig=1%, Sb=100, pb=0.1%, Change=5%", "Mig=1%, Sb=100, pb=0.1%, Change=25%"), ordered = TRUE)


options(scipen=10000)
dfe_plot = ggplot(expectation, aes(y=as.numeric(value), x=as.factor(coeff), shape=Population, color=Method)) +
  geom_point(position=position_dodge(width = .5),size=3) +
  # geom_errorbar(aes(ymin = CI5, ymax = CI95)) + 
  geom_pointrange(aes(ymin=CI5, ymax=CI95), position=position_dodge(width = .5)) +
  scale_color_brewer(palette="Dark2") +
  xlab("s") +
  ylab("Fraction") +
  #ylim(0,0.5) +
  #scale_y_log10(breaks=c(0.0001, 0.001, 0.01, 0.1, 0.5, 1)) +
  theme_bw() +  
  facet_wrap(~scenario, nrow = 6, ncol = 3) +
  theme(text = element_text(size=20), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom")
ggsave(plot = dfe_plot, "../plots/dadi/discretized_fullDFE.png", width=20, height=15, dpi="retina")



ben = subset(expectation, coeff == "> 0")
del = subset(expectation, coeff != "> 0")
ben$DFE = "Positive"
del$DFE = "Negative"
expectation = rbind(ben, del)

options(scipen=10000)
dfe_plot = ggplot(del, aes(y=as.numeric(value), x=as.factor(coeff), shape=Population, color=Method)) +
  geom_point(position=position_dodge(width = .5),size=3) +
  # geom_errorbar(aes(ymin = CI5, ymax = CI95)) + 
  geom_pointrange(aes(ymin=CI5, ymax=CI95), position=position_dodge(width = .5)) +
  scale_color_brewer(palette="Dark2") +
  xlab("s") +
  ylab("Fraction") +
  #ylim(0,0.5) +
  #scale_y_log10(breaks=c(0.0001, 0.001, 0.01, 0.1, 0.5, 1)) +
  theme_bw() +  
  facet_wrap(~scenario, nrow = 6, ncol = 3) +
  theme(text = element_text(size=20), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom", legend.title = element_blank())

ggsave(plot = dfe_plot, "../plots/dadi/discretized_delDFE.png", width=20, height=15, dpi="retina")

dfe_plot = ggplot(ben, aes(y=as.numeric(value), x=as.factor(coeff), shape=Population, color=Method)) +
  geom_point(position=position_dodge(width = .5),size=3) +
  # geom_errorbar(aes(ymin = CI5, ymax = CI95)) + 
  geom_pointrange(aes(ymin=CI5, ymax=CI95), position=position_dodge(width = .5)) +
  #scale_color_brewer(palette="Dark2") +
  xlab("s") +
  ylab("Fraction") +
  #ylim(0,0.5) +
  #scale_y_log10(breaks=c(0.0001, 0.001, 0.01, 0.1)) +
  theme_bw() +  
  facet_wrap(~scenario, scales = "free_y", nrow = 6, ncol = 3) +
  theme(text = element_text(size=20), axis.title.y=element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom")

ggsave(plot = dfe_plot, "../plots/dadi/discretized_benDFE.png", width=18, height=15, dpi="retina")



### Demographic parameters ###
args = commandArgs(trailingOnly = TRUE)

args[1] = 0
args[2] = 2
args[3] = 0.25

MIGRATION = args[1]
POSSEL = as.numeric(args[2])
CHANGE = args[3]

p <- vector('list', 18)
mig <- vector('list', 18)
pre <- vector('list', 18)
i=0
wb <- createWorkbook()
for (POSSEL in c(2, 20, 200))
{
  for (MIGRATION in c(0, 0.01))
  {
    for (CHANGE in c(0, 0.05, 0.25))
    {
      i = i+1
      print(paste(i, MIGRATION, POSSEL, CHANGE))
      dt = fread(paste0("../results/dadi_outputs/sim_", MIGRATION, "-", POSSEL, "-", CHANGE, "_demographic_confidence_intervals.txt") )
      colnames(dt) <- c("bound", "step_size", "Tpre", "nuPre", "Tdiv", "nu1div", "nu2div", "T1F", "T2F", "nu1F", "nu2F", "mw2d", "md2w", "missid")

      dt_step = data.frame(t(subset(dt, step_size == 1e-4)))
      #dt_step = subset(dt_step, X1 != "Lower_bound")
      
      dt_step$X1 = as.numeric(dt_step$X1)
      dt_step$X2 = as.numeric(dt_step$X2)

      dt_step$mean = ((dt_step$X1)+(dt_step$X2))/2
      #if (i == 13){
      #  dt_step$mean = c(NA, 0.0001, 0.073896373549883,	1.70449312822617,	0.511609212705555,	0.710288356088611,	0.505914407426143,	0.445655843912499,	0.421018882179885,	4.57471975300398,	4.65693251839453,	0.0147527753761547,	0.0109875137620283,	0.00247505568469394)
      #}
      
      max_TF = max(dt_step[c(8,9),]$mean)

      dt_step[dt_step < 0] <- 0

      wild_Ne_vector =   rev(c(dt_step[10,]$mean,dt_step[10,]$mean, dt_step[6,]$mean, dt_step[4,]$mean,1))
      dom_Ne_vector =    rev(c(dt_step[11,]$mean,dt_step[11,]$mean, dt_step[7,]$mean, dt_step[4,]$mean,1))

      wild_Ne_vector_low =   rev(c(dt_step[10,]$X1,dt_step[10,]$X1, dt_step[6,]$X1, dt_step[4,]$X1,1))
      wild_Ne_vector_up  =   rev(c(dt_step[10,]$X2,dt_step[10,]$X2, dt_step[6,]$X2, dt_step[4,]$X2,1))

      dom_Ne_vector_low  =    rev(c(dt_step[11,]$X1,dt_step[11,]$X1, dt_step[7,]$X1, dt_step[4,]$X1,1))
      dom_Ne_vector_up   =    rev(c(dt_step[11,]$X2,dt_step[11,]$X2, dt_step[7,]$X2, dt_step[4,]$X2,1))

      wild_Ne_true_vector = c(1,1)
      wild_time_true_vector = c(0,1)

      if (dt_step[8,]$mean == max_TF) {
        wild_time_vector = rev(c(0,dt_step[8,]$mean,  dt_step[8,]$mean+dt_step[5,]$mean,  dt_step[8,]$mean+dt_step[5,]$mean+dt_step[3,]$mean, (dt_step[8,]$mean+dt_step[5,]$mean+dt_step[3,]$mean)*1.5))
        wild_time_vector_low = rev(c(0,dt_step[8,]$X1,  dt_step[8,]$X1+dt_step[5,]$X1,  dt_step[8,]$X1+dt_step[5,]$X1+dt_step[3,]$X1, (dt_step[8,]$X1+dt_step[5,]$X1+dt_step[3,]$X1)*1.5))
        wild_time_vector_up  = rev(c(0,dt_step[8,]$X2,  dt_step[8,]$X2+dt_step[5,]$X2,  dt_step[8,]$X2+dt_step[5,]$X2+dt_step[3,]$X2, (dt_step[8,]$X2+dt_step[5,]$X2+dt_step[3,]$X2)*1.5))
        
        dom_time_vector =  rev(c(0,dt_step[9,]$mean,  dt_step[9,]$mean+dt_step[5,]$mean+(dt_step[8,]$mean-dt_step[9,]$mean), dt_step[9,]$mean+dt_step[5,]$mean+(dt_step[8,]$mean-dt_step[9,]$mean)+dt_step[3,]$mean,(dt_step[9,]$mean+dt_step[5,]$mean+(dt_step[8,]$mean-dt_step[9,]$mean)+dt_step[3,]$mean)*1.5))
        dom_time_vector_low =  rev(c(0,dt_step[9,]$X1,  dt_step[9,]$X1+dt_step[5,]$X1+(dt_step[8,]$X1-dt_step[9,]$X1), dt_step[9,]$X1+dt_step[5,]$X1+(dt_step[8,]$X1-dt_step[9,]$X1)+dt_step[3,]$X1,(dt_step[9,]$X1+dt_step[5,]$X1+(dt_step[8,]$X1-dt_step[9,]$X1)+dt_step[3,]$X1)*1.5))
        dom_time_vector_up  =  rev(c(0,dt_step[9,]$X2,  dt_step[9,]$X2+dt_step[5,]$X2+(dt_step[8,]$X2-dt_step[9,]$X2), dt_step[9,]$X2+dt_step[5,]$X2+(dt_step[8,]$X2-dt_step[9,]$X2)+dt_step[3,]$X2,(dt_step[9,]$X2+dt_step[5,]$X2+(dt_step[8,]$X2-dt_step[9,]$X2)+dt_step[3,]$X2)*1.5))
        
        } else {
        wild_time_vector = rev(c(0,dt_step[8,]$mean,  dt_step[8,]$mean+dt_step[5,]$mean+(dt_step[9,]$mean-dt_step[8,]$mean), dt_step[8,]$mean+dt_step[5,]$mean+(dt_step[9,]$mean-dt_step[8,]$mean)+dt_step[3,]$mean,(dt_step[8,]$mean+dt_step[5,]$mean+(dt_step[9,]$mean-dt_step[8,]$mean)+dt_step[3,]$mean)*1.5))
        wild_time_vector_low = rev(c(0,dt_step[8,]$X1,  dt_step[8,]$X1+dt_step[5,]$X1+(dt_step[9,]$X1-dt_step[8,]$X1), dt_step[8,]$X1+dt_step[5,]$X1+(dt_step[9,]$X1-dt_step[8,]$X1)+dt_step[3,]$X1,(dt_step[8,]$X1+dt_step[5,]$X1+(dt_step[9,]$X1-dt_step[8,]$X1)+dt_step[3,]$X1)*1.5))
        wild_time_vector_up  = rev(c(0,dt_step[8,]$X2,  dt_step[8,]$X2+dt_step[5,]$X2+(dt_step[9,]$X2-dt_step[8,]$X2), dt_step[8,]$X2+dt_step[5,]$X2+(dt_step[9,]$X2-dt_step[8,]$X2)+dt_step[3,]$X2, (dt_step[8,]$X2+dt_step[5,]$X2+(dt_step[9,]$X2-dt_step[8,]$X2)+dt_step[3,]$X2)*1.5))
        
        dom_time_vector  = rev(c(0,dt_step[9,]$mean,  dt_step[9,]$mean+dt_step[5,]$mean,  dt_step[9,]$mean+dt_step[5,]$mean+dt_step[3,]$mean, ( dt_step[9,]$mean+dt_step[5,]$mean+dt_step[3,]$mean)*1.5))
        dom_time_vector_low  = rev(c(0,dt_step[9,]$X1,  dt_step[9,]$X1+dt_step[5,]$X1,  dt_step[9,]$X1+dt_step[5,]$X1+dt_step[3,]$X1, (dt_step[9,]$X1+dt_step[5,]$X1+dt_step[3,]$X1)*1.5))
        dom_time_vector_up   = rev(c(0,dt_step[9,]$X2,  dt_step[9,]$X2+dt_step[5,]$X2,  dt_step[9,]$X2+dt_step[5,]$X2+dt_step[3,]$X2, (dt_step[9,]$X2+dt_step[5,]$X2+dt_step[3,]$X2)*1.5))
      }

      yLIM = max(c(dom_Ne_vector,wild_Ne_vector))+0.1
      xLIM = max(c(dom_time_vector,wild_time_vector))

      dom_Ne_true_vector =   rev(c(1,   1,   0.04,1   ,1))
      dom_time_true_vector = rev(c(xLIM,0.20,0.20,0.18,0))
      
      wild = cbind(wild_time_vector, wild_Ne_vector, "Wild", "Inferred")
      dom = cbind(dom_time_vector, dom_Ne_vector, "Domestic", "Inferred")
      dom2 = cbind(dom_time_true_vector, dom_Ne_true_vector, "Domestic", "True")
      dt_plot = data.frame(rbind(wild, dom, dom2))
      colnames(dt_plot) <- c("Time", "Ne", "Population", "State")
      dt_plot$Time = as.numeric(dt_plot$Time)
      dt_plot$Ne = as.numeric(dt_plot$Ne)
      
      p[[i]] <- ggplot(dt_plot, aes(x=Time, y=Ne, color=Population, linetype=State)) +
        geom_step(direction = "vh", size=2, alpha=0.6) +
        scale_color_brewer(palette="Dark2") +
        theme_bw() +
        theme(legend.position="none", plot.title = element_text(hjust = 0.5, size = 8)) +
        ggtitle(paste0("ID = ", i, "\nm = ", MIGRATION, ", Sb = ", POSSEL/2, " and pc = ", CHANGE))

      #pdf(paste0("../plots/Demography.pdf"), width = 10, height = 6) 
      #pdf(paste0("../plots/Demography_", MIGRATION, "-", POSSEL, "-", CHANGE, ".pdf"), width = 4, height = 3) 
      #par(mar=c(1,1,1,1))
      #par( mfrow= c(6,3) ) 
      #plot(wild_time_vector, wild_Ne_vector, ylim = c(0,yLIM), xlim = c(0,xLIM), type='s', col="#994F00", 
      #     xlab = "Time", ylab = "Ne", lwd = 3, main = (paste0("Migration = ", MIGRATION, ", Sb = ", POSSEL/2, " and pc = ", CHANGE)))
      #par(new=TRUE)
      #plot(dom_time_vector, dom_Ne_vector,   ylim = c(0,yLIM), xlim = c(0,xLIM), type="s", col="#006CD1", xlab = "", ylab = "", lwd = 3, xaxt="n", yaxt="n")
      #par(new=TRUE)
      #plot(dom_time_true_vector, dom_Ne_true_vector,   ylim = c(0,yLIM), xlim = c(0,xLIM), type="s", col="black", xlab = "Time", ylab = "Ne", lwd = 2, lty=3)
      #dev.off() 
      #plot(dom_time_vector_low, dom_Ne_vector_low,   ylim = c(0,yLIM), xlim = c(0,xLIM), type="s", col="#006CD1", xlab = "Time", ylab = "Ne", lwd = 2, lty=3)
      #par(new=TRUE)
      #plot(dom_time_vector_up, dom_Ne_vector_up,     ylim = c(0,yLIM), xlim = c(0,xLIM), type="s", col="#006CD1", xlab = "", ylab = "", lwd = 2, lty=3, xaxt="n", yaxt="n")
      #par(new=TRUE)
      #plot(wild_time_vector_low, wild_Ne_vector_low, ylim = c(0,yLIM), xlim = c(0,xLIM), type="s", col="#994F00", xlab = "", ylab = "", lwd = 2, lty=3, xaxt="n", yaxt="n")
      #par(new=TRUE)
      #plot(wild_time_vector_up, wild_Ne_vector_up,   ylim = c(0,yLIM), xlim = c(0,xLIM), type="s", col="#994F00", xlab = "", ylab = "", lwd = 2, lty=3, xaxt="n", yaxt="n")
      #par(new=TRUE)
      #plot(dom_time_true_vector, dom_Ne_true_vector,   ylim = c(0,yLIM), xlim = c(0,xLIM), type="s", col="black", xlab = "", ylab = "", lwd = 2, lty=3, xaxt="n", yaxt="n")
    
      dt_step = dt_step[-c(1,2),]
      colnames(dt_step) <- c("CI_2.5%", "CI_97.5%", "value")
      
      sheet_name = paste0("SCENARIO_ID=", i)
      addWorksheet(wb,sheet_name)
      writeData(wb, sheet_name, dt_step, rowNames=TRUE)
      
      dt_step$params = rownames(dt_step)

      ### Pre Ne plots ###
      dt_step_pre <- subset(dt_step, params == "Tpre" | params == "nuPre") %>% 
        mutate(epoch = case_when(`CI_2.5%` == 0 & params == "Tpre" ~  "Not Significant",
                                (`CI_2.5%` <= 1 & `CI_97.5%` >= 1) & params == "nuPre" ~  "Not Significant"))
      dt_step_pre["epoch"][is.na(dt_step_pre["epoch"])] <- "Significant"

      pre[[i]] <- ggplot(dt_step_pre, aes(x=params, y=value, color=epoch)) +
        geom_point() +
        geom_errorbar(aes(ymin=`CI_2.5%`, ymax=`CI_97.5%`), width=.2,
                      position=position_dodge(0.05)) +
        #ylim(0,2) +
        xlab("") + ylab("") +
        facet_wrap(~params, scales = "free") +
        scale_color_manual(breaks = c("Significant", "Not Significant"),
                           values=c("red", "grey")) +
        theme_bw() +
        theme(legend.position="none", plot.title = element_text(hjust = 0.5, size = 8)) +
        ggtitle(paste0("ID = ", i, "\nm = ", MIGRATION, ", Sb = ", POSSEL/2, " and pc = ", CHANGE))
      
      
      ### Migration plots ###
      dt_step_mig <- subset(dt_step, params == "md2w" | params == "mw2d") %>% 
        mutate(epoch = case_when(`CI_2.5%` == 0 & params == "md2w" ~  "Not Significant",
                                 `CI_2.5%` == 0 & params == "mw2d" ~  "Not Significant"))
      dt_step_mig["epoch"][is.na(dt_step_mig["epoch"])] <- "Significant"
      
      mig[[i]] <- ggplot(dt_step_mig, aes(x=params, y=value, color=epoch)) +
        geom_point() +
        geom_errorbar(aes(ymin=`CI_2.5%`, ymax=`CI_97.5%`), width=.2,
                      position=position_dodge(0.05)) +
        #ylim(0,2) +
        xlab("") + ylab("") +
        #facet_wrap(~params) +
        scale_color_manual(breaks = c("Significant", "Not Significant"),
                           values=c("red", "grey")) +
        theme_bw() +
        geom_hline(yintercept = 0.01) +
        geom_hline(yintercept = 0, color = "white") +
        
        theme(legend.position="none", plot.title = element_text(hjust = 0.5, size = 8)) +
        ggtitle(paste0("ID = ", i, "\nm = ", MIGRATION, ", Sb = ", POSSEL/2, " and pc = ", CHANGE))
    }
  }
}

saveWorkbook(wb, "../results/Supp_Table1.xlsx", overwrite = TRUE)

#pervasive_and_nearlyneutral = grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], nrow = 2)
#common_and_weak = grid.arrange(p[[7]], p[[8]], p[[9]], p[[10]], p[[11]], p[[12]], nrow = 2)
#rare_and_strong = grid.arrange(p[[13]], p[[14]], p[[15]], p[[16]], p[[17]], p[[18]], nrow = 2)

all = grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], p[[7]], p[[8]], p[[9]], p[[10]], p[[11]], p[[12]], p[[13]], p[[14]], p[[15]], p[[16]], p[[17]], p[[18]], nrow = 6)
ggsave(all, filename = "../plots/All_Scenarios_Demography.png",width = 8.5, height = 12, dpi = 300)

all = grid.arrange(pre[[1]], pre[[2]], pre[[3]], pre[[4]], pre[[5]], pre[[6]], pre[[7]], pre[[8]], pre[[9]], pre[[10]], pre[[11]], pre[[12]], pre[[13]], pre[[14]], pre[[15]], pre[[16]], pre[[17]], pre[[18]], nrow = 6)
ggsave(all, filename = "../plots/All_Scenarios_nuPre.png",width = 8.5, height = 12, dpi = 300)

all = grid.arrange(mig[[1]], mig[[2]], mig[[3]], mig[[4]], mig[[5]], mig[[6]], mig[[7]], mig[[8]], mig[[9]], mig[[10]], mig[[11]], mig[[12]], mig[[13]], mig[[14]], mig[[15]], mig[[16]], mig[[17]], mig[[18]], nrow = 6)
ggsave(all, filename = "../plots/All_Scenarios_Migration.png",width = 8.5, height = 12, dpi = 300)

