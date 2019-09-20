                   
setwd('F:/Delia/COMPARE/dtangle.data_0.3.6.tar/dtangle.data/data')
source('https://raw.githubusercontent.com/jumphone/Delia/master/Delia.R')










##################
load('abbas.rda')
load('becht.rda')
load('gong.rda')
load('kuhn.rda')
load('linsley.rda')
load('liu.rda')
load('newman_fl.rda')
load('newman_pbmc.rda')
load('shen_orr.rda')
load('shi.rda')
load('siegert.rda')




DATA.ALL=list()

DATA.ALL[[1]]=abbas
DATA.ALL[[2]]=becht
DATA.ALL[[3]]=gong
DATA.ALL[[4]]=kuhn
DATA.ALL[[5]]=linsley
DATA.ALL[[6]]=liu
DATA.ALL[[7]]=newman_fl
DATA.ALL[[8]]=newman_pbmc
DATA.ALL[[9]]=shen_orr
DATA.ALL[[10]]=shi
DATA.ALL[[11]]=siegert



DELIA.COR=c()
DTAN.COR=c()
TOP=2000


data.index=1
while(data.index<=length(DATA.ALL)){
    DATA=DATA.ALL[[data.index]]

  
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

DELIA.COR=c(DELIA.COR, mean(COR))

##############################################


DT.RESULT=dtangle(Y=DATA$data$log, reference=t(REF),n_markers=TOP )
PP.DT=DT.RESULT$estimates[MIX.INDEX,]

COR=c()
i=1
while(i<=ncol(REF)){
    this_cor=cor(PP.TRUE[,i],PP.DT[,i])
    COR=c(COR,this_cor)
    i=i+1
   }

DTAN.COR=c(DTAN.COR, mean(COR))




#################

  
 print(data.index)
  
    data.index=data.index+1

   }

































