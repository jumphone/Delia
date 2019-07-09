# Second Benchmark:


# PBMC data derived from https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html

source('Delia.R')
library(Seurat)


pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(pbmc,'./RDS/pbmcdemo.RDS')

PBMC.REF=.generate_ref(pbmc@assays$RNA@data, as.character(pbmc@active.ident))
saveRDS(PBMC.REF,'./RDS/PBMC.REF.RDS')
###############################


REF=PBMC.REF
REF=apply(REF, 2, .norm_exp)


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
saveRDS(ALLR, file='./RDS/PBMC.ALLR.RDS')
saveRDS(EXP, file='./RDS/PBMC.EXP.RDS')
##############################


# Non-linear transfer 
getRanMap<-function(x){
  oldx=x
  x[order(x)]=sort(rnorm(length(x)))
  y=x-min(x)
  y[which(oldx==0)]=0
  return(y)
  }

set.seed(1)
REXP=apply(EXP,2,getRanMap)


# Add Random Noise
addNOI=function(x){
     M=mean(x)
     y=x+M/3 * rnorm(length(x))
     return(y)
  }

REXP=REXP
set.seed(123)
REXP=t(apply(REXP,1,addNOI))
REXP[which( REXP<0)]=0
rownames( REXP)=rownames(REXP)
colnames(REXP)=colnames(REXP)

REXP=apply(REXP, 2, .norm_exp)
colnames(REXP)=paste0('NOI_',c(1:ncol(REXP)))
colnames(EXP)=paste0('RAN_',c(1:ncol(EXP)))
##############################
saveRDS(REXP, file='./RDS/PBMC.REXP.RDS')
##############################




mydelia=Delia(REXP, REF ,COMBAT=TRUE, WEIGHT=TRUE)
# 4s
RATIO=mydelia$out
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')

CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,3]+CORMAT[4,4]+CORMAT[5,5]+
CORMAT[6,6]+CORMAT[7,7]+CORMAT[8,8]+CORMAT[9,9]
#ACS = 7.48

library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))


#######################################



# MuSiC_all, 32 seconds
library(MuSiC)
library('Biobase')
library(xbioc)
library(pvar)
###########################
sc_mat=pbmc@assays$RNA@data #read.table('Zeisel_exp_sc_mat.txt',sep='\t',header=T,row.names=1)
sc_mat=apply(sc_mat, 2, .norm_exp)
TAG= as.character(pbmc@active.ident) #read.table('Zeisel_exp_sc_mat_cluster_merged.txt',header=T,sep='\t')
CELLTYPE=TAG
assayData=REXP #readRDS('./RDS/REXP.RDS')

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



























