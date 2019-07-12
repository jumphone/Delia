<img src="https://github.com/jumphone/Delia/blob/master/img/Delia_LOGO.png" width="200">


### An ultrafast cell type deconvolution tool

DELIA: (D) (E) convo (L) ut (I) on (A) pproach

Author: Feng Zhang

Environment: R 

</br>

# Requirement:

Please install "ComBat":

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
        
    BiocManager::install("sva")
    BiocManager::install("limma")

</br>

# Usage:

* [1. Load package](#1-Load-package)
* [2. Load DEMO data](#2-Load-DEMO-data)
* [3. Use single-cell data as the reference](#3-Use-single-cell-data-as-the-reference)
* [4. Run Delia](#4-Run-Delia)
* [5. Results](#5-Results)
* [6. Visualization](#6-Visualization)

</br>


## 1. Load package:

    library('sva')
    library('limma')
    
Directly load online file:

    source('https://raw.githubusercontent.com/jumphone/Delia/master/Delia.R')
    
or, download "Delia.R" and load it:
    
    source('Delia.R')
    
</br>

## 2. Load DEMO data:

Download DEMO data from: https://github.com/jumphone/Delia/tree/master/DEMO

    sc_exp_mat <- .readTable(PATH='SingleCellREF/sc_exp_mat.txt', SEP='\t')
    sc_exp_mat=log(sc_exp_mat+1,10)
    sc_exp_mat=apply(sc_exp_mat, 2, .norm_exp)
    
    tag <- read.table('SingleCellREF/tag.txt', sep='\t', header=TRUE,row.names=1)
    tag <- as.character(tag[,1])
    
    EXP <- readRDS('EXP.RDS')  # This matrix has been log-normalized already. 
    true_ratio <- readRDS('TRUE_RATIO.RDS')

    
</br>

## 3. Use single-cell data as the reference:
    
    REF <- .generate_ref(sc_exp_mat, tag)
    
sc_exp_mat: single-cell expression matrix (row: gene, col: cell)

tag: cell type labels (a character vector of cell types)

</br>

## 4. Run Delia:

#### 4.1 Linear Regression (Manuscript Version):

    mydelia <- Delia(EXP, REF, COMBAT=TRUE)      

EXP: expression matrix of query data; colname is query name, rowname is gene name

REF: expression matrix of reference; colname is cell type, rowname is gene name

COMBAT: use ComBat to do batch-effect correction (default is TRUE)

</br>

#### 4.2 Principal Components Regression (PCR):

#### If there are some very similar cell types (Multicollinearity) in your reference,  please try PCR

#### Caution: PCR will be slower than linear regression.

    install.packages('pls')    
    install.packages('tcltk2')

    mydelia <- Delia(EXP, REF, COMBAT=TRUE, PCR=TRUE, PCV=0.95, SHOW=TRUE)   

PCR: use Principal Components Regression (PCR) (default is FALSE). 

PCV: when PCR is TRUE, set cutoff for variance explained by used PCs (default is 0.95)

SHOW: show progress bar (default is FALSE). 

</br>

## 5. Results:   

#### 5.1 Proportion Matrix:

If cell types in the query and reference are <b> WELL MATCHED </b>, please use "mydelia$out".

For each query sample, sum is normalized to 1.

Colname is query name; Rowname is cell type.

    mydelia$out

</br>

#### 5.2 Coefficient Matrix:

If you're <b> NOT CLEAR ABOUT THE CELL TYPES </b> in your query data, please use "mydelia$coef".
 
For each query sample, sum is not normalized to 1.

Colname is query name; Rowname is cell type.

    mydelia$coef
    
</br>

## 6. Visualization:   
   
    estimate_ratio <- mydelia$out    
    estimate_ratio_coef <- mydelia$coef
    
    # Please use "mydelia$coef" when cell types in query and reference are not well matched.
    
    # Save output (TXT format):   
    .writeTable(DATA=estimate_ratio, PATH='OUTPUT.txt', SEP='\t')
    .writeTable(DATA=estimate_ratio_coef, PATH='OUTPUT_COEF.txt', SEP='\t')


#### 6.1 Show results of the first 10 query samples:

#### Proportion Plot

    show_ratio <- estimate_ratio[,1:10]
    
    library('gplots')
     
    heatmap.2(t(show_ratio),scale=c("none"), dendrogram='none',
        Rowv=F,Colv=F,cellnote=round(t(show_ratio),2), notecol='black',
        trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),
        margins=c(10,10))


<img src="https://raw.githubusercontent.com/jumphone/Delia/master/img/PLOT1.png" width="300">
  

    COLOR=colorRampPalette(c('grey90','grey30'))(nrow(show_ratio))
   
    show_ratio_bar=cbind(show_ratio,rep(0,nrow(show_ratio)) )
   
    barplot(show_ratio_bar,las=2,col=COLOR)
    legend("topright", 
           legend = rownames(show_ratio), 
           fill = COLOR)
       
<img src="https://raw.githubusercontent.com/jumphone/Delia/master/img/PLOT1_BAR.png" width="350">
   
   
#### Coefficient Plot

    show_ratio_coef <- estimate_ratio_coef[,1:10]
    
    library('gplots')
     
    heatmap.2(t(show_ratio_coef),scale=c("none"), dendrogram='none',
        Rowv=F,Colv=F,cellnote=round(t(show_ratio_coef),2), notecol='black',
        trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),
        margins=c(10,10))
  
<img src="https://raw.githubusercontent.com/jumphone/Delia/master/img/PLOT1_COEF.png" width="300">
    


#### 6.2 Show correlation between estimated and true ratios of "A":

#### Proportion Plot

    True.A=true_ratio[1,]
    Est.A=estimate_ratio[1,]
    PCC=round(cor(True.A,Est.A,method='pearson'),2)
    
    plot(True.A, Est.A, pch=16, main=paste0('PCC=',PCC))
           
<img src="https://raw.githubusercontent.com/jumphone/Delia/master/img/PLOT2.png" width="300">
        

#### Coefficient Plot

    True.A=true_ratio[1,]
    Est.A=estimate_ratio_coef[1,]
    PCC=round(cor(True.A,Est.A,method='pearson'),2)
    
    plot(True.A, Est.A, pch=16, main=paste0('PCC=',PCC))
           
<img src="https://raw.githubusercontent.com/jumphone/Delia/master/img/PLOT2_COEF.png" width="300">



