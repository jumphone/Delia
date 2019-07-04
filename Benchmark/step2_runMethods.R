# Load Data
REF=readRDS('./RDS/REF.RDS')
ALLR=readRDS('./RDS/ALLR.RDS')
EXP=readRDS('./RDS/EXP.RDS')
REXP=readRDS('./RDS/REXP.RDS')
SC.REF=readRDS('./RDS/SC.REF.RDS')
V.SC.REF=readRDS('./RDS/V.SC.REF.RDS')
V.SC.REF_scmat=readRDS('./RDS/V.SC.REF_scmat.RDS')
################

# Load Package 
library(sva)
library(limma)
source('Delia.R')
################


# Delia with all genes, 2 seconds
mydelia = Delia(REXP, SC.REF, COMBAT=TRUE) 
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
#########################


# Delia with variable genes, less than 1 second
mydelia = Delia(REXP, V.SC.REF, COMBAT=TRUE) 
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


# CIBERSORT with variable genes, 144 seconds
CB=read.table('./RESULT/CIBERSORT.Output_Job3.txt',header=T,row.names=1,sep='\t')
RATIO=t(CB[,c(1:(ncol(CB)-3))])
#######
saveRDS(RATIO, file='./RESULT/CIBERSORT_var.RDS')
#######
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
pdf('./RESULT/CIBERSORT_var.pdf',width=6, height=6)
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
dev.off()
#########################


# CIBERSORTx with variable genes, 274 seconds
CB=read.table('./RESULT/CIBERSORTx_Job2_Adjusted.txt',header=T,row.names=1,sep='\t')
RATIO=t(CB[,c(1:(ncol(CB)-3))])
#######
saveRDS(RATIO, file='./RESULT/CIBERSORTx_var.RDS')
#######
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
pdf('./RESULT/CIBERSORTx_var.pdf',width=6, height=6)
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
dev.off()
#########################








