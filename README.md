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

#### Step1. Load pacakge:

    library('sva')
    library('limma')
    source('https://raw.githubusercontent.com/jumphone/Delia/master/Delia.R')
    # source('Delia.R')
    
#### Step2. If you are using single-cell data as reference, please do the following step:
    
    REF=.generate_ref(sc_exp_mat, tag)
    
sc_exp_mat: single-cell expression matrix (row: gene, col: cell)

tag: cell type labels (a character vector of cell types)

#### Step3. Run Delia:
    
    mydelia=Delia(EXP, REF, COMBAT=TRUE)
 
#### Input:            

EXP：expression matrix of query data; colname is query name, rowname is gene

REF: expression matrix of reference; colname is cell type, rowname is gene 

COMBAT: use ComBat to do batch-effect correction (default is TRUE)

#### Result:   

mydelia$out: sum is normalized to 1 (proportion matrix); colname is query name, rowname is cell type

mydelia$coef: sum is not normalized (coefficient matrix); colname is query name, rowname is cell type
    
If cell types in the query and reference are well matched, please use "mydelia$out"   

If you're not clear about the cell types in your query data, please use "mydelia$coef"
    



