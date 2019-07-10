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

## 2. Use single-cell data as the reference:
    
    REF=.generate_ref(sc_exp_mat, tag)
    
sc_exp_mat: single-cell expression matrix (row: gene, col: cell)

tag: cell type labels (a character vector of cell types)

</br>

## 3. Run Delia:
    
    mydelia=Delia(EXP, REF, COMBAT=TRUE)      

EXP: expression matrix of query data; colname is query name, rowname is gene

REF: expression matrix of reference; colname is cell type, rowname is gene 

COMBAT: use ComBat to do batch-effect correction (default is TRUE)

</br>

## 4. Results:   

#### 4.1 Proportion Matrix

If cell types in the query and reference are well matched, please use "mydelia$out".

For each query sample, sum is normalized to 1.

Colname is query name, and rowname is cell type.

    mydelia$out

#### 4.2 Coefficient Matrix

If you're not clear about the cell types in your query data, please use "mydelia$coef".
 
For each query sample, sum is not normalized to 1.

Colname is query name, and rowname is cell type.

    mydelia$coef
    
    
## 5. Visualization:   
   

   



