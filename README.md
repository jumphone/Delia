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
    
    # sc_exp_mat: log-normalized (non-scale) single-cell expression matrix (row: gene, col: cell)
    # tag: cell type labels
    
    
    mydelia=Delia(EXP, REF, COMBAT=TRUE)
        
    # Input:            
    # EXP：colname is sample name, rowname is gene
    # REF: colname is cell type, rowname is gene 
    # The expression values of EXP and REF should be log normalized   
    # COMBAT: use ComBat to do batch-effect correction (default is TRUE)

    # Result:   
    # mydelia$out: colname is sample name, rowname is cell type
    
    



