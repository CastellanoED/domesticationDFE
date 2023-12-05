args = commandArgs(trailingOnly=TRUE)
if(length(args)<3) {
  args[1] <- 1000 #3000 #10000 #1000 #5000
  args[2] <- 150 #600 #1500 #150 #600 
  args[3] <- "../results/fastas/REF_A.fa" #~/Desktop"
}
nloci <- as.numeric(args[1])
size_gene <- as.numeric(args[2])
fileRefA <- args[3]

len <- nloci * size_gene
seqs <- array("A",dim=c(1,len))
#WRITE FASTA FILE
write.table(x=sprintf(">REFERENCE"),file=sprintf("%s",fileRefA),append=F,quote=F,sep="",eol="\n",row.names=F,col.names=F)
write.table(x=seqs,file=sprintf("%s",fileRefA),append=T, quote=F,sep="",eol="",row.names=F,col.names=F)
write.table(x="\n",file=sprintf("%s",fileRefA),append=T, quote=F,sep="",eol="",row.names=F,col.names=F)

show("run finished")
