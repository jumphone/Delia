# Brain Bulk Data, used to generate simulated bulk data
# Ref: Zhang, Y., et al. An RNA-sequencing transcriptome and splicing database of glia, neurons, and vascular cells of the cerebral cortex. J Neurosci 2014;34(36):11929-11947.

REF=read.table('Reference_expression.txt',sep='\t',header=T,row.names=1)

REF[,4] = REF[,4] + REF[,5]
REF=REF[,c(1,2,3,4,6,7)]

REF=log(REF+1,10)
REF=apply(REF, 2, .norm_exp)
hist(REF[,1])

colnames(REF)=c('ASTRO','NEURON','OPC','OLIGO','MICRO','ENDO')



###############
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
###############

EXP=apply(EXP, 2, .norm_exp)


#####################
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
