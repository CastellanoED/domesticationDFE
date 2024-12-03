library("data.table")
library("dplyr")
library("tidyr")
library("ggplot2")
library(gridExtra)
library(openxlsx)


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

dt_best_loglk$ancestral_Ne <- dt_best_loglk$theta_nsyn / 80
write.table(dt_best_loglk, file = "../results/data_frames/dadi_best_fits.txt", row.names = F, col.names = T, sep = "\t", quote = F)





# tapply(dt_best_loglk$ancestral_Ne, dt_best_loglk$params, mean)
# table(list(dt_best_loglk$true_possel, dt_best_loglk$assumed_possel, dt_best_loglk$true_pchange, dt_best_loglk$mig))
# tapply(dt_best_loglk$pchange, dt_best_loglk$params, quantile, probs=c(0.025,0.975))

# # mu*L is 20 for nsyn
# (100 * 2/3  * 120   * 10000    * 2.5e-7)  * 4 * 5000 #expected theta nsyn under neutrality
# #(it * nsyn * #loci * bp/locus * mu/site)  * heritable units * "2"  * N
# #The "2" comes from the fact that two sequences that have diverged for time t are different by 2 * mu * t mutations, since both diverging lineages accumulate mutations.

# ### Comparison of Ne inference across populations and  methods ###
# demog = fread("../results/data_frames/flexible_dem_best_fits.txt", fill = T)
# colnames(demog) <- c("params", "loglk", "Tpre", "nuPre", "Tdiv", "nu1div", "nu2div", "T1F", "T2F", "nu1F", "nu2F", "mw2d", "md2w", "missid","theta_syn", "ratio_nsyn_syn")
# demog$Wild = (demog$nu1F*demog$theta_syn) / 40
# demog$Domesticated = (demog$nu2F*demog$theta_syn) / 40
# demog$Ancestral = (demog$theta_syn) / 40
# demog_Ne = dplyr::select(demog, "params", "Ancestral", "Wild", "Domesticated")

# tapply(demog_Ne$Ancestral, demog_Ne$params, mean)
# demog_Ne <- gather(demog_Ne, Population, Ne, Ancestral:Domesticated, factor_key=TRUE)
# demog_Ne$Method = "dadi"
# demog_Ne_polyDFE = fread("../results/data_frames/polyDFE_Ne.txt", header = T)
# demog_Ne = rbind(demog_Ne, demog_Ne_polyDFE)


