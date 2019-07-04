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
mydelia = Delia(REXP, SC.REF, COMBAT=TRUE, WEIGHT=TRUE) 
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
mydelia = Delia(REXP, V.SC.REF, COMBAT=TRUE, WEIGHT=TRUE) 
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





RATIO=readRDS(file='./RESULT/Delia_all.RDS')
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
#4.16

RATIO=readRDS(file='./RESULT/Delia_var.RDS')
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
#3.90


RATIO=readRDS(file='./RESULT/CIBERSORT_var.RDS')
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
#3.49

RATIO=readRDS(file='./RESULT/CIBERSORTx_var.RDS')
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
#3.67

RATIO=readRDS(file='./RESULT/MuSiC_all.RDS')
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
#4.02


RATIO=readRDS(file='./RESULT/MuSiC_var.RDS')
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]
#3.43








