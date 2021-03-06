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
    
<b>sc_exp_mat</b>: single-cell expression matrix (row: gene, col: cell)

<b>tag</b>: cell type labels (a character vector of cell types)

</br>

## 4. Run Delia:

    mydelia <- Delia(EXP, REF, COMBAT=TRUE, METHOD='lm', SHOW=FALSE, RANK=FALSE)      

<b>EXP</b>: expression matrix of query data; Colname is query name, rowname is gene name

<b>REF</b>: expression matrix of reference; Colname is cell type, rowname is gene name

<b>COMBAT</b>: use ComBat to do batch-effect correction (Default is TRUE)

<b>SHOW</b>: show progress bar ('tcltk2' package) (Default is FALSE)

<b>RANK</b>: use Rank-based normalization (Default if FASLE)
    
<b>METHOD</b>: 
    
* <b>'lm'</b> &nbsp; : &nbsp; Linear model (used in our manuscript) (Default).
    
* <b>'rlm'</b> &nbsp; : &nbsp; Robust linear model ('MASS' package). 
    
* <b>'pcr'</b> &nbsp; : &nbsp; Principal components regression ('pls' package). Caution: 'pcr' is much slower than 'lm' and 'rlm'.

* <b>'opt'</b> &nbsp; : &nbsp; Linear constraint optimization ('optim' function). RANK will be changed to TRUE when using 'opt'.
    
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

#### 5.3 Save Result:   

RDS format:

    saveRDS(mydelia, file='mydelia.RDS')
    
TXT format:

    .writeTable(DATA=mydelia$out, PATH='out.txt', SEP='\t')
    .writeTable(DATA=mydelia$coef, PATH='coef.txt', SEP='\t')

</br>


## 6. Visualization:   


#### 6.1 Show results of the first 10 query samples:

#### Proportion:

    show_ratio <-  mydelia$out[,1:10]
    
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
   
   
#### Coefficient:

    show_ratio_coef <- mydelia$coef[,1:10]
    
    library('gplots')
     
    heatmap.2(t(show_ratio_coef),scale=c("none"), dendrogram='none',
        Rowv=F,Colv=F,cellnote=round(t(show_ratio_coef),2), notecol='black',
        trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),
        margins=c(10,10))
  
<img src="https://raw.githubusercontent.com/jumphone/Delia/master/img/PLOT1_COEF.png" width="300">
 
Normalize coefficients:

     show_ratio_coef <- .scaleGlobal(mydelia$coef)[,1:10]
     
     library('gplots')
     
     heatmap.2(t(show_ratio_coef),scale=c("none"), dendrogram='none',
        Rowv=F,Colv=F,cellnote=round(t(show_ratio_coef),2), notecol='black',
        trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),
        margins=c(10,10))

<img src="https://raw.githubusercontent.com/jumphone/Delia/master/img/PLOT1_COEF_SCALE.png" width="300">
 
 
#### 6.2 Show correlation between estimated and true ratios of "A":

#### Proportion:

    True.A=true_ratio[1,]
    Est.A=mydelia$out[1,]
    PCC=round(cor(True.A,Est.A,method='pearson'),2)
    
    plot(True.A, Est.A, pch=16, main=paste0('PCC=',PCC))
           
<img src="https://raw.githubusercontent.com/jumphone/Delia/master/img/PLOT2.png" width="300">
        

#### Coefficient:

    True.A=true_ratio[1,]
    Est.A=mydelia$coef[1,]
    PCC=round(cor(True.A,Est.A,method='pearson'),2)
    
    plot(True.A, Est.A, pch=16, main=paste0('PCC=',PCC))
           
<img src="https://raw.githubusercontent.com/jumphone/Delia/master/img/PLOT2_COEF.png" width="300">



