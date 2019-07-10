<img src="https://github.com/jumphone/Delia/blob/master/img/Delia_LOGO.png" width="200">


### An ultrafast cell type deconvolution tool

DELIA: (D) (E) convo (L) ut (I) on (A) pproach

Author: Feng Zhang

Environment: R 

</br>

# Requirement:

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
        
    BiocManager::install("sva")
    BiocManager::install("limma")

</br>

# Usage:

## 1. Load pacakge:

    library('sva')
    library('limma')
    
Directly load online file:

    source('https://raw.githubusercontent.com/jumphone/Delia/master/Delia.R')
    
or, download "Delia.R" and load it:
    
    source('Delia.R')
    
</br>

## 2. Load DEMO data:

Download DEMO data from: https://github.com/jumphone/Delia/tree/master/DEMO

    sc_exp_mat <- read.table('SingleCellREF/sc_exp_mat.txt', sep='\t',row.names=1,header=T)
    sc_exp_mat=log(sc_exp_mat+1,10)
    sc_exp_mat=apply(sc_exp_mat, 2, .norm_exp)
    
    tag <- read.table('SingleCellREF/tag.txt', sep='\t',row.names=1,header=T)
    tag <- as.character(tag[,1])
    
    EXP <- readRDS('EXP.RDS')    
    true_ratio <- readRDS('TRUE_RATIO.RDS')

    
</br>

## 3. Use single-cell data as the reference:
    
    REF <- .generate_ref(sc_exp_mat, tag)
    
sc_exp_mat: single-cell expression matrix (row: gene, col: cell)

tag: cell type labels (a character vector of cell types)

</br>

## 4. Run Delia:
    
    mydelia <- Delia(EXP, REF, COMBAT=TRUE)      

EXP: expression matrix of query data; colname is query name, rowname is gene

REF: expression matrix of reference; colname is cell type, rowname is gene 

COMBAT: use ComBat to do batch-effect correction (default is TRUE)

</br>

## 5. Results:   

#### 5.1 Proportion Matrix

If cell types in the query and reference are well matched, please use "mydelia$out".

For each query sample, sum is normalized to 1.

Colname is query name, and rowname is cell type.

    mydelia$out

#### 5.2 Coefficient Matrix

If you're not clear about the cell types in your query data, please use "mydelia$coef".
 
For each query sample, sum is not normalized to 1.

Colname is query name, and rowname is cell type.

    mydelia$coef
    
</br>

## 5. Visualization:   
   
    estimate_ratio <- mydelia$out 
    
    
Show first 10 demo query samples:

    show_rato <- estimate_ratio[,1:10]
    
    library('gplots')
     
    heatmap.2(t(show_rato),scale=c("none"), dendrogram='none',
        Rowv=F,Colv=F,cellnote=round(t(show_rato),2), notecol='black',
        trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),
        margins=c(10,10))




Show correlation between estimated and true ratios of "A":
    
    True.A=true_ratio[1,]
    Est.A=estimate_ratio[1,]
    PCC=round(cor(True.A,Est.A,method='pearson'),2)
    plot(True.A, Est.A, pch=16, main=paste0('PCC=',PCC))
    
        
        
       


