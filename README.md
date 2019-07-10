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

# Usage:

    library('sva')
    library('limma')
    source('https://raw.githubusercontent.com/jumphone/Delia/master/Delia.R')
    # source('Delia.R')
    
    # If you are using single-cell data as reference, please do the following step:
    
    REF=.generate_ref(sc_exp_mat, tag)
    
    # sc_exp_mat: single-cell expression matrix (row: gene, col: cell)
    # tag: cell type labels
    
    
    mydelia=Delia(EXP, REF, COMBAT=TRUE)
        
    # Input:            
    # EXP：colname is query name, rowname is gene
    # REF: colname is cell type, rowname is gene 
    # COMBAT: use ComBat to do batch-effect correction (default is TRUE)

    # Result:   
    # mydelia$out: sum is normalized to 1 (proportion matrix); colname is query name, rowname is cell type
    # mydelia$coef: sum is not normalized (coefficient matrix); colname is query name, rowname is cell type
    
    # If you're not clear about the cell types in your query data, please use "mydelia$coef"
    



