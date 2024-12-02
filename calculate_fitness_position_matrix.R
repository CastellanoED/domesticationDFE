args = commandArgs(trailingOnly=TRUE)

if(length(args)<10) {
  args[1] <- 1828105181250
  args[2] <- 10000
  args[3] <- 120
  args[4] <- 0.002
  args[5] <- -0.02
  args[6] <- 0.3
  args[7] <- 0.75
  args[8] <- 0.25
  args[9] <- 0.75
  args[10] <- c("../results/fitness_position_matrix/fitness_position_matrix.txt")
} 
seed.number <- as.numeric(args[1])
nloci <- as.numeric(args[2])
size_gene <- as.numeric(args[3])
s_mean_beneficial <- as.numeric(args[4])
s_mean_deleterious <- as.numeric(args[5])
shape_deleterious <- as.numeric(args[6])
prop_del_anc <- as.numeric(args[7])
change_prop <- as.numeric(args[8])
prop_del_new <- as.numeric(args[9])
file_fm <- (args[10])

set.seed(seed.number)

pos_deleterious_no_modified <- prop_del_anc * (1.0 - change_prop);
pos_deleterious_changed_del <- prop_del_anc * change_prop * prop_del_new;
pos_deleterious_changed_ben <- prop_del_anc * change_prop * (1.0 - prop_del_new);

pos_beneficial_no_modified <- (1.0 - prop_del_anc) * (1.0 - change_prop);
pos_beneficial_changed_del <- (1.0 - prop_del_anc) * change_prop * prop_del_new;
pos_beneficial_changed_ben <- (1.0 - prop_del_anc) * change_prop * (1.0 - prop_del_new);

fm=matrix(rep(0.0,5*nloci*size_gene),nrow=nloci*size_gene,ncol=5,byrow=F);
colnames(fm) <- c("Pos","Type","s.p1","s.p2","s.p3")
vec <- rep(c(1,1,0),nloci*size_gene/3)
tm  <- sample(c(2:7),size=nloci*size_gene,replace=T,
             prob=c(pos_beneficial_no_modified,pos_beneficial_changed_ben,pos_beneficial_changed_del,
                    pos_deleterious_no_modified,pos_deleterious_changed_del,pos_deleterious_changed_ben));
tm <- tm * vec
tm[(tm==0)] <- 1
#Position
fm[,1] <- 1:(nloci*size_gene)
#type of mutations
fm[,2] <- tm
#ancestral beneficial
fm[fm[,2]>=2 & fm[,2]<=4,3] <- rexp(sum(fm[,2]>=2 & fm[,2]<=4),1/s_mean_beneficial)
fm[fm[,2]>=2 & fm[,2]<=4,4] <- fm[fm[,2]>=2 & fm[,2]<=4,3]
fm[fm[,2]>=2 & fm[,2]<=4,5] <- fm[fm[,2]>=2 & fm[,2]<=4,4]
fm[fm[,2]==3,5] <- rexp(sum(fm[,2]==3),1/s_mean_beneficial)
fm[fm[,2]==4,5] <- -rgamma(sum(fm[,2]==4),scale=-s_mean_deleterious/shape_deleterious,shape=shape_deleterious)
#ancestral deleterious
fm[fm[,2]>=5 & fm[,2]<=7,3] <- -rgamma(sum(fm[,2]>=5 & fm[,2]<=7),scale=-s_mean_deleterious/shape_deleterious,shape=shape_deleterious)
fm[fm[,2]>=5 & fm[,2]<=7,4] <- fm[fm[,2]>=5 & fm[,2]<=7,3]
fm[fm[,2]>=5 & fm[,2]<=7,5] <- fm[fm[,2]>=5 & fm[,2]<=7,4]
fm[fm[,2]==6,5] <- -rgamma(sum(fm[,2]==6),scale=-s_mean_deleterious/shape_deleterious,shape=shape_deleterious)
fm[fm[,2]==7,5] <- rexp(sum(fm[,2]==7),1/s_mean_beneficial)

#keep matrix into a file
write.table(x=fm,file=file_fm,quote=F,col.names=T,row.names=F,append=F,sep="\t")
show("run finished")
