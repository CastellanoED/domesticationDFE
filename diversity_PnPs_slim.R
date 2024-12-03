library("data.table")
library("dplyr")
library("tidyr")
library("ggplot2")
library(gridExtra)
library(openxlsx)
#library(ggpubr)


### Calculating diversity and Pn/Ps across independent simulation run ###
pNpS_all_sims = data.frame()

for (it in seq(1, 100))
{
  for (mig in c(0, 0.01))
  {
    for (possel in c(0, 2, 20, 200))
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
pNpS_all_sims$params = paste0("m = ", pNpS_all_sims$mig, ",", "Sb = ", pNpS_all_sims$possel, " and ", "pc = ", pNpS_all_sims$pchange)

A <- ggplot(pNpS_all_sims, aes(x=ratio, fill=pop ) ) + 
  geom_histogram(alpha=0.5, position="identity") +
  facet_wrap(~ params, ncol = 3) +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  xlab("Pn/Ps") +
  theme(text = element_text(size=15), axis.title.y=element_blank(),legend.position="none")
A

armonic_number = sum(1/seq(39))
L = 120/3 * 10000

B <- ggplot(pNpS_all_sims, aes(x=pS/L/armonic_number, fill=pop ) ) + 
  geom_histogram(alpha=0.5, position="identity") +
  facet_wrap(~ params, ncol = 3) +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  xlab(expression(theta[syn])) +
  theme(text = element_text(size=15), axis.title.y=element_blank(),legend.position="none")
B

ggsave(B, file="../plots/Supp_Figure6.png", width=12, height=15, dpi=200)
ggsave(A, file="../plots/Supp_Figure7.png", width=12, height=15, dpi=200)
