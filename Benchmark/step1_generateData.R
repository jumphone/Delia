source('Delia.R')

###################################################################
# Mouse Brain Bulk Data, used to generate simulated bulk data
# https://github.com/jumphone/scRef/tree/master/Reference/MouseBrain_Bulk_Zhang2014
# Ref: Zhang, Y., et al. An RNA-sequencing transcriptome and splicing database of glia, neurons, and vascular cells of the cerebral cortex. J Neurosci 2014;34(36):11929-11947.
#######################################################################

# Load Data
REF=read.table('Reference_expression.txt',sep='\t',header=T,row.names=1)
REF[,4] = REF[,4] + REF[,5]
REF=REF[,c(1,2,3,4,6,7)]
REF=log(REF+1,10)
REF=apply(REF, 2, .norm_exp)
#hist(REF[,1])
colnames(REF)=c('ASTRO','NEURON','OPC','OLIGO','MICRO','ENDO')
##############################
saveRDS(REF, file='./RDS/REF.RDS')
##############################

# Generate Simulated Data
NUM=100
ALLR=matrix(0,ncol=NUM,nrow=ncol(REF))
EXP=matrix(0,ncol=NUM,nrow=nrow(REF))

set.seed(1)
i=1
while(i<=NUM){
    this_r=rnorm(ncol(REF))
    this_r=pnorm(this_r)
  
    #this_r=runif(this_r)
    this_r=this_r/sum(this_r)
    this_exp=REF %*% this_r
    ALLR[,i]=this_r
    EXP[,i]=this_exp
    if(i%%10==1){print(i)}
i=i+1
}

rownames(ALLR)=colnames(REF)
rownames(EXP)=rownames(REF)
EXP=apply(EXP, 2, .norm_exp)
##############################
saveRDS(ALLR, file='./RDS/ALLR.RDS')
saveRDS(EXP, file='./RDS/EXP.RDS')
##############################


# Add Random Noise
getRanMap<-function(x){
  oldx=x
  x[order(x)]=sort(rnorm(length(x)))
  y=x-min(x)
  y[which(oldx==0)]=0
  return(y)
  }

set.seed(1)
REXP=apply(EXP,2,getRanMap)

REXP=REXP
  set.seed(123)
  addNOI=function(x){
     M=mean(x)
     y=x+M/3 * rnorm(length(x))#(runif(length(x))*2-1)
     return(y)
  }

  REXP=t(apply(REXP,1,addNOI))
  REXP[which( REXP<0)]=0
  rownames( REXP)=rownames(REXP)
  colnames(REXP)=colnames(REXP)

REXP=apply(REXP, 2, .norm_exp)
colnames(REXP)=paste0('NOI_',c(1:ncol(REXP)))
colnames(EXP)=paste0('RAN_',c(1:ncol(EXP)))
##############################
saveRDS(REXP, file='./RDS/REXP.RDS')
##############################



###################################################################
# Mouse Brain single-cell Data, used to generate reference
# https://github.com/jumphone/scRef/tree/master/Benchmark
# Ref: Zeisel et al. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq, Science
#######################################################################

# Load Data
sc_mat=read.table('Zeisel_exp_sc_mat.txt',sep='\t',header=T,row.names=1)
sc_mat=log(sc_mat+1,10)
sc_mat=apply(sc_mat, 2, .norm_exp)
TAG=read.table('Zeisel_exp_sc_mat_cluster_merged.txt',header=T,sep='\t')
table(TAG[,2])
TAG[,2]=as.character(TAG[,2])
TAG[which(TAG[,2]=='astrocytes_ependymal'),2]='ASTRO'
TAG[which(TAG[,2]=='endothelial-mural'),2]='ENDO'
TAG[which(TAG[,2]=='microglia'),2]='MICRO'
TAG[which(TAG[,2]=='neurons'),2]='NEURON'
TAG[which(TAG[,2]=='oligodendrocytes'),2]='OLIGO'

# Reference with All Gene ##############
SC.REF=.generate_ref(sc_mat, TAG, min_cell=1)
SC.REF=apply(SC.REF, 2, .norm_exp)
SC.REF=SC.REF[,c(1,4,5,3,2)]
##############################
saveRDS(SC.REF, file='./RDS/SC.REF.RDS')
##############################

# Reference with Variable Gene##############
VAR=apply(SC.REF, 1, var)
VG=which(rank(-VAR) <= 2000  )
V.SC.REF=SC.REF[VG,]
##############################
saveRDS(V.SC.REF, file='V.SC.REF.RDS')
##############################

# Single-cell Reference Matrix with Varible Gene##############
V.SC.REF_scmat=sc_mat[VG,]
colnames(V.SC.REF_scmat)=TAG[,2]
##############################
saveRDS(V.SC.REF_scmat, file='./RDS/V.SC.REF_scmat.RDS')
##############################



# Generate CIBERSORT files
##################
OUT=cbind(rownames(EXP),EXP)
colnames(OUT)[1]='GENE'
write.table(OUT, file='./CIBERSORT/EXP_mix.txt',sep='\t',row.names=F,col.names=T,quote=F)
##################

##################
OUT=cbind(rownames(REXP),REXP)
colnames(OUT)[1]='GENE'
write.table(OUT, file='./CIBERSORT/REXP_mix.txt',sep='\t',row.names=F,col.names=T,quote=F)
##################

##################
OUT=cbind(rownames(V.SC.REF),V.SC.REF)
colnames(OUT)[1]='GENE'
write.table(OUT, file='./CIBERSORT/V.SC.REF_sig.txt',sep='\t',row.names=F,col.names=T,quote=F)
##################

##################
OUT=cbind(rownames(V.SC.REF_scmat),V.SC.REF_scmat)
colnames(OUT)[1]='GENE'
write.table(OUT, file='./CIBERSORT/V.SC.REF_scmat_sig.txt',sep='\t',row.names=F,col.names=T,quote=F)
##################

