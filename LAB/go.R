setwd('F:/Delia/LAB')
source('Delia.R')

library(Seurat)
pbmc.data <- Read10X(data.dir = "./filtered_gene_matrix/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features =VariableFeatures(object = pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)

VEC=pbmc@reductions$umap@cell.embeddings
K=kmeans(VEC,centers=50)
C=K$cluster
Idents(pbmc)=as.character(C)

#pbmc <- FindNeighbors(pbmc, dims = 1:10)
#pbmc <- FindClusters(pbmc, resolution = 0.5)

DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "pbmc.RDS")


FeaturePlot(pbmc,features=c('CCL3','TMEM119'))
##########################################################



#########################################################
pbmc=readRDS('pbmc.RDS')


#EXP=.generate_ref(as.matrix(pbmc@assays$RNA@counts), Idents(pbmc))
EXP=.generate_ref(as.matrix(pbmc@assays$RNA@data), Idents(pbmc))
#EXP=log(EXP+1,10)


setwd('F:/Delia/LAB')
source('Delia.R')
library(Seurat)

REF <- .readTable(PATH='CMAT.txt.rmCC.txt', SEP='\t', UP=TRUE)
REF=log(REF+1,10)

########################
source('Delia.R')
mydelia <- Delia(EXP, REF, COMBAT=TRUE)   



#################

show_ratio_coef <- mydelia$coef

library('gplots') 
heatmap.2(t(show_ratio_coef),scale=c("none"), dendrogram='none',
    Rowv=F,Colv=F,cellnote=round(t(show_ratio_coef),2), notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),
    margins=c(10,10))


#################




##################
show_ratio <- mydelia$out

library('gplots')
heatmap.2(t(show_ratio),scale=c("none"), dendrogram='none',
    Rowv=F,Colv=F,cellnote=round(t(show_ratio),2), notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),
    margins=c(10,10))



show_ratio=apply(show_ratio,2,scale)

library('gplots')
heatmap.2(t(show_ratio),scale=c("none"), dendrogram='none',
    Rowv=F,Colv=F,cellnote=round(t(show_ratio),2), notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),
    margins=c(10,10))




