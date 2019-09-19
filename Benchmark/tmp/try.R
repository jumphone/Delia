                   
setwd('F:/Delia/COMPARE/dtangle.data_0.3.6.tar/dtangle.data/data')
source('https://raw.githubusercontent.com/jumphone/Delia/master/Delia.R')





##################
load('abbas.rda')
DATA=abbas
##############################################

EXP=t(DATA$data$log)

REF.FLAG=apply(DATA$annotation$pure,1,sum)
which(REF.FLAG==1)

REF.TAG=rep(NA,length(REF.FLAG))
i=1
while(i<=nrow(DATA$annotation$pure)){
    if(REF.FLAG[i]==1){
        this_tag=colnames(DATA$annotation$pure)[which(DATA$annotation$pure[i,]==TRUE)]
        REF.TAG[i]=this_tag
        }

    i=i+1}

REF.INDEX=which(REF.FLAG==1)

REF=.generate_ref(EXP[,REF.INDEX], REF.TAG[REF.INDEX])


MIX.INDEX=which(REF.FLAG==0)
MIX=EXP[,MIX.INDEX]

##############

library(dtangle)

MK=find_markers(Y=DATA$data$log, reference=t(REF) )

TOP=1000

MK.GENE=c()
i=1
while(i<=ncol(REF)){
    MK.GENE=c(MK.GENE, names(MK$V[[i]])[1:TOP])
i=i+1}

###################
REF.MK=REF[which(rownames(REF) %in% MK.GENE),]

mydelia <- Delia(MIX, REF.MK, COMBAT=FALSE, METHOD='lm', SHOW=FALSE, RANK=FALSE)     
PP.TRUE=DATA$annotation$mixture[MIX.INDEX,]
PP.DELIA=t(mydelia$out)

PP.TRUE=PP.TRUE[,order(colnames(PP.TRUE))]
PP.DELIA=PP.DELIA[,order(colnames(PP.DELIA))]

COR=c()
i=1
while(i<=ncol(REF)){
    this_cor=cor(PP.TRUE[,i],PP.DELIA[,i])
    COR=c(COR,this_cor)
    i=i+1
   }
mean(COR)
#0.9670429

##############################################

TOP=1000
DT.RESULT=dtangle(Y=DATA$data$log, reference=t(REF),n_markers=TOP )
PP.DT=DT.RESULT$estimates[MIX.INDEX,]

COR=c()
i=1
while(i<=ncol(REF)){
    this_cor=cor(PP.TRUE[,i],PP.DT[,i])
    COR=c(COR,this_cor)
    i=i+1
   }
mean(COR)
#0.9617639



#################






##################
load('becht.rda')
DATA=becht
##############################################

EXP=t(DATA$data$log)

REF.FLAG=apply(DATA$annotation$pure,1,sum)
which(REF.FLAG==1)

REF.TAG=rep(NA,length(REF.FLAG))
i=1
while(i<=nrow(DATA$annotation$pure)){
    if(REF.FLAG[i]==1){
        this_tag=colnames(DATA$annotation$pure)[which(DATA$annotation$pure[i,]==TRUE)]
        REF.TAG[i]=this_tag
        }

    i=i+1}

REF.INDEX=which(REF.FLAG==1)

REF=.generate_ref(EXP[,REF.INDEX], REF.TAG[REF.INDEX])


MIX.INDEX=which(REF.FLAG==0)
MIX=EXP[,MIX.INDEX]

##############

library(dtangle)

MK=find_markers(Y=DATA$data$log, reference=t(REF) )

TOP=1000

MK.GENE=c()
i=1
while(i<=ncol(REF)){
    MK.GENE=c(MK.GENE, names(MK$V[[i]])[1:TOP])
i=i+1}

###################
REF.MK=REF[which(rownames(REF) %in% MK.GENE),]

mydelia <- Delia(MIX, REF.MK, COMBAT=TRUE, METHOD='lm', SHOW=FALSE, RANK=FALSE)     
PP.TRUE=DATA$annotation$mixture[MIX.INDEX,]
PP.DELIA=t(mydelia$out)

PP.TRUE=PP.TRUE[,order(colnames(PP.TRUE))]
PP.DELIA=PP.DELIA[,order(colnames(PP.DELIA))]

COR=c()
i=1
while(i<=ncol(REF)){
    this_cor=cor(PP.TRUE[,i],PP.DELIA[,i])
    COR=c(COR,this_cor)
    i=i+1
   }
mean(COR)
#0.7711682

##############################################

TOP=1000
DT.RESULT=dtangle(Y=DATA$data$log, reference=t(REF),n_markers=TOP )
PP.DT=DT.RESULT$estimates[MIX.INDEX,]

COR=c()
i=1
while(i<=ncol(REF)){
    this_cor=cor(PP.TRUE[,i],PP.DT[,i])
    COR=c(COR,this_cor)
    i=i+1
   }
mean(COR)
#0.7421834



#################










