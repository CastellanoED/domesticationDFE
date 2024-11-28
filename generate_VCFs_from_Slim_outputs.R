library("data.table")
library("dplyr")
library("tidyr")

args = commandArgs(trailingOnly=TRUE)

MIGRATION = args[1]
POSSEL = args[2]
CHANGE = args[3]
ITER = args[4]

wild = fread(paste0("../results/Simulations/sim_", MIGRATION, "-", POSSEL, "-", CHANGE, "-", ITER, "_slim_output_sample_file.txt_pop1"))
dom  = fread(paste0("../results/Simulations/sim_", MIGRATION, "-", POSSEL, "-", CHANGE, "-", ITER, "_slim_output_sample_file.txt_pop2"))

wild = dplyr::select(wild, V2, V3, V4, V9)
colnames(wild) <- c("id", "type", "pos", "freq")
wild$pop = "wild"

dom = dplyr::select(dom, V2, V3, V4, V9)
colnames(dom) <- c("id", "type", "pos", "freq")
dom$pop = "dom"

genomes = t(data.frame(numeric(40)))

dom = cbind(dom, genomes)
dom_long <- gather(dom, genome_ID, genotype, V1:V40, factor_key=TRUE)
dom_long$genome_ID = gsub("V", "D", dom_long$genome_ID)

wild = cbind(wild, genomes)
wild_long <- gather(wild, genome_ID, genotype, V1:V40, factor_key=TRUE)
wild_long$genome_ID = gsub("V", "W", wild_long$genome_ID)

both = rbind(wild_long, dom_long)

ids = unique(both$id)
positions = unique(both$pos)

recurrent_mutations = length(ids) - length(positions)
print(paste("Number of recurrent mutations =", recurrent_mutations))

vcf_like = data.frame()
for (ID in ids)
{
  # ID = "147485659"
  site = subset(both, id == ID)
  wild_site = subset(site, pop == "wild")
  dom_site  = subset(site, pop == "dom")
  
  if (nrow(wild_site)==0)
  {
    wild_site = data.frame(id=ID, type=dom_site$type, pos=dom_site$pos, freq=0, pop="wild", genome_ID=paste0("W",seq(40)), genotype=0)
  }
  
  if (nrow(dom_site)==0)
  {
    dom_site = data.frame(id=ID, type=wild_site$type, pos=wild_site$pos, freq=0, pop="dom", genome_ID=paste0("D",seq(40)), genotype=0)
  }
  
  wild_site_sample = sample_n(wild_site, size=unique(wild_site$freq), replace = FALSE)
  if (nrow(wild_site_sample)>0) {wild_site_sample$genotype = 1}
  wild_site_antisample = anti_join(wild_site, wild_site_sample, by = c("genome_ID"))
  wild_site = rbind(wild_site_sample, wild_site_antisample)
  
  dom_site_sample = sample_n(dom_site, size=unique(dom_site$freq), replace = FALSE)
  if (nrow(dom_site_sample)>0) {dom_site_sample$genotype = 1}
  dom_site_antisample = anti_join(dom_site, dom_site_sample, by = c("genome_ID"))
  dom_site = rbind(dom_site_sample, dom_site_antisample)
  
  wild_site = dplyr::select(wild_site, -freq, -pop)
  dom_site = dplyr::select(dom_site, -freq, -pop)
  
  both_genot = rbind(dom_site, wild_site)
  both_genot_wide <- spread(both_genot, genome_ID, genotype)
  vcf_like = rbind(vcf_like, both_genot_wide)
  
}

sorted_vcf_like = vcf_like[order(vcf_like$pos),]
# #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT

`#CHROM` = ITER
POS = sorted_vcf_like$pos
ID = sorted_vcf_like$id
REF = "A"
ALT = "T"
QUAL = "."
FILTER = "."
INFO = "AA=A"
FORMAT = "GT"

vcf = cbind(`#CHROM`, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sorted_vcf_like)

VCF_syn  = subset(vcf, type=="m1")
VCF_nsyn = subset(vcf, type=="m2")

VCF_syn  = dplyr::select(VCF_syn,  -type, -pos, -id)
VCF_nsyn = dplyr::select(VCF_nsyn, -type, -pos, -id)

write.table(VCF_syn,  paste0("../results/VCF/sim_", MIGRATION, "-", POSSEL, "-", CHANGE, "-", ITER, "_syn.vcf"),  sep="\t", append=F, col.names = T, row.names = F, quote = F)
write.table(VCF_nsyn, paste0("../results/VCF/sim_", MIGRATION, "-", POSSEL, "-", CHANGE, "-", ITER, "_nsyn.vcf"), sep="\t", append=F, col.names = T, row.names = F, quote = F)
