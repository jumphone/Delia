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
#4.163812

##########################################
#2019.9.26
library(dtangle)
TOP=2000
CCC=.simple_combine(REXP,SC.REF)
DT.RESULT=dtangle(Y=t(CCC$exp_sc_mat1), reference=t(CCC$exp_sc_mat2),n_markers=TOP )

PP.DT=t(DT.RESULT$estimates)

CORMAT.DT=cor(t(PP.DT), t(ALLR), method='pearson')


pdf('./RESULT/dtangle_top2000.pdf',width=6, height=6)
library('gplots')
heatmap.2(CORMAT.DT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT.DT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
dev.off()


CORMAT.DT[1,1]+CORMAT.DT[2,2]+CORMAT.DT[3,4]+CORMAT.DT[3,5]+CORMAT.DT[4,6]+CORMAT.DT[5,7]
#3.361094

##################################

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




# 20190712###########
source('Delia.R')
mydelia = Delia(REXP, SC.REF)
RATIO=mydelia$out
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
######





# CIBERSORT with variable genes, 145 seconds
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


# CIBERSORTx with variable genes, 275 seconds
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




# MuSiC_all, 32 seconds
library(MuSiC)
library('Biobase')
library(xbioc)
library(pvar)
###########################
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
CELLTYPE=TAG[,2]
assayData=readRDS('./RDS/REXP.RDS')

Sys.time()
Bulk.eset=ExpressionSet(assayData)
sampleID=1:ncol(sc_mat)
SubjectName=colnames(sc_mat)
cellTypeID=CELLTYPE
cellType=CELLTYPE
pheno.matrix=data.frame(sampleID,SubjectName,cellTypeID,cellType)
rownames(pheno.matrix)=colnames(sc_mat)
metadata <- data.frame(labelDescription= c("Sample ID", "Subject Name", "Cell Type ID", "Cell Type Name"), row.names=c("sampleID", "SubjectName", "cellTypeID", "cellType"))
SC.eset = ExpressionSet(assayData = data.matrix(sc_mat), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix, varMetadata = metadata) )
sampleType=colnames(sc_mat)
PP_ALL = music_prop(bulk.eset = Bulk.eset, sc.eset = SC.eset ,markers = NULL, clusters = cellType, samples = sampleType, verbose = F)
Sys.time()

RATIO=t(PP_ALL$Est.prop.weighted)
RATIO=RATIO[c(5,1,2,3,4),]
#######
saveRDS(RATIO, file='./RESULT/MuSiC_all.RDS')
#######
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
pdf('./RESULT/MuSiC_all.pdf',width=6, height=6)
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
dev.off()
#########################



# MuSiC_var, 6 seconds
library(MuSiC)
library('Biobase')
library(xbioc)
library(pvar)
###########################
V.SC.REF_scmat=readRDS('./RDS/V.SC.REF_scmat.RDS')
assayData=readRDS('./RDS/REXP.RDS')
colnames(V.SC.REF_scmat)=colnames(sc_mat)
sc_mat=V.SC.REF_scmat

Sys.time()
Bulk.eset=ExpressionSet(assayData)
sampleID=1:ncol(sc_mat)
SubjectName=colnames(sc_mat)
cellTypeID=CELLTYPE
cellType=CELLTYPE
pheno.matrix=data.frame(sampleID,SubjectName,cellTypeID,cellType)
rownames(pheno.matrix)=colnames(sc_mat)
metadata <- data.frame(labelDescription= c("Sample ID", "Subject Name", "Cell Type ID", "Cell Type Name"), row.names=c("sampleID", "SubjectName", "cellTypeID", "cellType"))
SC.eset = ExpressionSet(assayData = data.matrix(sc_mat), phenoData =  new("AnnotatedDataFrame", data = pheno.matrix, varMetadata = metadata) )
sampleType=colnames(sc_mat)
PP_ALL = music_prop(bulk.eset = Bulk.eset, sc.eset = SC.eset ,markers = NULL, clusters = cellType, samples = sampleType, verbose = F)
Sys.time()


RATIO=t(PP_ALL$Est.prop.weighted)
RATIO=RATIO[c(5,1,2,3,4),]
#######
saveRDS(RATIO, file='./RESULT/MuSiC_var.RDS')
#######
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
pdf('./RESULT/MuSiC_var.pdf',width=6, height=6)
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
dev.off()
#########################











ALLR=readRDS('./RDS/ALLR.RDS')

#install.packages('PerformanceAnalytics')
library(PerformanceAnalytics)

###################


RATIO=readRDS(file='./RESULT/Delia_all.RDS')
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
#4.16

CALLR=cbind(ALLR[1,], ALLR[2,], ALLR[4,]+ALLR[5,], ALLR[6,],ALLR[7,])
TRATIO=t(RATIO)
colnames(CALLR)=colnames(TRATIO)
colnames(CALLR)[3]='NFOL.add.MOL'
MAT=cor(CALLR, TRATIO)
MAT[1,1]+MAT[2,2]+MAT[3,3]+MAT[4,4]+MAT[5,5]
#4.03

colnames(CALLR)=paste0('True.',colnames(CALLR))
colnames(TRATIO)=paste0('Estimate.',colnames(TRATIO))
DATA=cbind(CALLR, TRATIO)
tiff('./RESULT/Delia_cor_all.tiff',width=800,height=800)
chart.Correlation(DATA, histogram=TRUE, pch=16)
dev.off()
#################################


RATIO=readRDS(file='./RESULT/Delia_var.RDS')
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
#3.90

CALLR=cbind(ALLR[1,], ALLR[2,], ALLR[4,]+ALLR[5,], ALLR[6,],ALLR[7,])
TRATIO=t(RATIO)
colnames(CALLR)=colnames(TRATIO)
colnames(CALLR)[3]='NFOL.add.MOL'
MAT=cor(CALLR, TRATIO)
MAT[1,1]+MAT[2,2]+MAT[3,3]+MAT[4,4]+MAT[5,5]
#3.76

colnames(CALLR)=paste0('True.',colnames(CALLR))
colnames(TRATIO)=paste0('Estimate.',colnames(TRATIO))
DATA=cbind(CALLR, TRATIO)
tiff('./RESULT/Delia_cor_var.tiff',width=800,height=800)
chart.Correlation(DATA, histogram=TRUE, pch=16)
dev.off()

############################################

RATIO=readRDS(file='./RESULT/CIBERSORT_var.RDS')
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
#3.49

CALLR=cbind(ALLR[1,], ALLR[2,], ALLR[4,]+ALLR[5,], ALLR[6,],ALLR[7,])
TRATIO=t(RATIO)
colnames(CALLR)=colnames(TRATIO)
colnames(CALLR)[3]='NFOL.add.MOL'
MAT=cor(CALLR, TRATIO)
MAT[1,1]+MAT[2,2]+MAT[3,3]+MAT[4,4]+MAT[5,5]
#3.36

colnames(CALLR)=paste0('True.',colnames(CALLR))
colnames(TRATIO)=paste0('Estimate.',colnames(TRATIO))
DATA=cbind(CALLR, TRATIO)
tiff('./RESULT/CIBERSORT_cor_var.tiff',width=800,height=800)
chart.Correlation(DATA, histogram=TRUE, pch=16)
dev.off()

########################################

RATIO=readRDS(file='./RESULT/CIBERSORTx_var.RDS')
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
#3.67

CALLR=cbind(ALLR[1,], ALLR[2,], ALLR[4,]+ALLR[5,], ALLR[6,],ALLR[7,])
TRATIO=t(RATIO)
colnames(CALLR)=colnames(TRATIO)
colnames(CALLR)[3]='NFOL.add.MOL'
MAT=cor(CALLR, TRATIO)
MAT[1,1]+MAT[2,2]+MAT[3,3]+MAT[4,4]+MAT[5,5]
#3.53

colnames(CALLR)=paste0('True.',colnames(CALLR))
colnames(TRATIO)=paste0('Estimate.',colnames(TRATIO))
DATA=cbind(CALLR, TRATIO)
tiff('./RESULT/CIBERSORTx_cor_var.tiff',width=800,height=800)
chart.Correlation(DATA, histogram=TRUE, pch=16)
dev.off()


########################################


RATIO=readRDS(file='./RESULT/MuSiC_all.RDS')
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
#4.02

CALLR=cbind(ALLR[1,], ALLR[2,], ALLR[4,]+ALLR[5,], ALLR[6,],ALLR[7,])
TRATIO=t(RATIO)
colnames(CALLR)=colnames(TRATIO)
colnames(CALLR)[3]='NFOL.add.MOL'
MAT=cor(CALLR, TRATIO)
MAT[1,1]+MAT[2,2]+MAT[3,3]+MAT[4,4]+MAT[5,5]
#3.87

colnames(CALLR)=paste0('True.',colnames(CALLR))
colnames(TRATIO)=paste0('Estimate.',colnames(TRATIO))
DATA=cbind(CALLR, TRATIO)
tiff('./RESULT/MuSiC_cor_all.tiff',width=800,height=800)
chart.Correlation(DATA, histogram=TRUE, pch=16)
dev.off()



##########################################

RATIO=readRDS(file='./RESULT/MuSiC_var.RDS')
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
#3.43

CALLR=cbind(ALLR[1,], ALLR[2,], ALLR[4,]+ALLR[5,], ALLR[6,],ALLR[7,])
TRATIO=t(RATIO)
colnames(CALLR)=colnames(TRATIO)
colnames(CALLR)[3]='NFOL.add.MOL'
MAT=cor(CALLR, TRATIO)
MAT[1,1]+MAT[2,2]+MAT[3,3]+MAT[4,4]+MAT[5,5]
#3.30

colnames(CALLR)=paste0('True.',colnames(CALLR))
colnames(TRATIO)=paste0('Estimate.',colnames(TRATIO))
DATA=cbind(CALLR, TRATIO)
tiff('./RESULT/MuSiC_cor_var.tiff',width=800,height=800)
chart.Correlation(DATA, histogram=TRUE, pch=16)
dev.off()

#######################################






