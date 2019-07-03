<img src="https://github.com/jumphone/Delia/blob/master/img/Delia_LOGO.png" width="200">


### Delia: an ultrafast cell-type deconvelution approach

Author: Feng Zhang

Environment: R 

</br>

# Requirement

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
        
    BiocManager::install("sva")
    BiocManager::install("limma")

# Usage


    source('https://raw.githubusercontent.com/jumphone/Delia/master/Delia.R')
    
    mydelia=Delia(EXP, REF, COMBAT=TRUE)
        
    # EXP：colname is sample name, rowname is gene
    # REF: colname is cell type, rowname is gene
    # COMBAT: use ComBat to do batch-effect correction (default is TRUE)

    # mydelia$out: colname is sample name, rowname is cell type
    
    



