setwd('F:/Delia')
# Load Data
REF=readRDS('./RDS/REF.RDS')
ALLR=readRDS('./RDS/ALLR.RDS')
EXP=readRDS('./RDS/EXP.RDS')
REXP=readRDS('./RDS/REXP.RDS')
SC.REF=readRDS('./RDS/SC.REF.RDS')
V.SC.REF=readRDS('./RDS/V.SC.REF.RDS')
V.SC.REF_scmat=readRDS('./RDS/V.SC.REF_scmat.RDS')
################
VAR=apply(SC.REF, 1, var)
VG=which(rank(-VAR) <= 2000  )
########################

library(deconvSeq)
############################


sc_mat=read.table('Zeisel_exp_sc_mat.txt',sep='\t',header=T,row.names=1)
sc_mat[which(is.na(sc_mat))]
SCMAT=sc_mat[VG,]

colnames(SCMAT)=paste0(colnames(SCMAT),'_',c(1:ncol(SCMAT)))
tissue.sc = colnames(V.SC.REF_scmat)
design.sc = model.matrix(~-1+as.factor(tissue.sc))
colnames(design.sc) = levels(as.factor(tissue.sc))
rownames(design.sc) = names(tissue.sc)
dge.sc = getdge(SCMAT, design.sc, ncpm.min=1, nsamp.min=4, method="bin.loess")
b0.sc = getb0.rnaseq(dge.sc, design.sc, ncpm.min=1, nsamp.min=4)
##########################################


source('Delia.R')
COM=.simple_combine(REXP, SCMAT)$com
BATCH=c(rep('EXP',ncol(REXP)), rep('REF',ncol(SCMAT)))
COM.combat=.combat(COM, BATCH) 


####################################################
dge.tissue = getdge(REXP,design=NULL, ncpm.min=1,nsamp.min=4)
####################################################
resultx1_s.sc = getx1.rnaseq(NB0=50, b0.sc, dge.tissue)






# Load Package 
library(sva)
library(limma)
################

source('Delia.R')
# Delia with all genes, 5 seconds
mydelia = Delia(REXP, SC.REF) 
RATIO=mydelia$out
#######
saveRDS(RATIO, file='./RESULT/Delia_all.RDS')
#######
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
pdf('./RESULT/Delia_all.pdf',width=6, height=6)
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
dev.off()
########################



RATIO=readRDS(file='./RESULT/Delia_all.RDS')
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]

plot(RATIO[4,],ALLR[6,],pch=16)


# Delia with variable genes, 1 second
mydelia = Delia(REXP, V.SC.REF) 
RATIO=mydelia$out
#######
saveRDS(RATIO, file='./RESULT/Delia_var.RDS')
#######
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
pdf('./RESULT/Delia_var.pdf',width=6, height=6)
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
dev.off()
#########################

RATIO=readRDS(file='./RESULT/Delia_var.RDS')
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]








