<img src="https://github.com/jumphone/Delia/blob/master/img/Delia_LOGO.png" width="200">


Delia: DEconvoLutIon Approach

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
    
    mydelia=Delia(EXP, REF, COMBAT=TRUE, WEIGHT=TRUE)
        
    # The expression values of EXP and REF should be log normalized   
       
    # EXP：colname is sample name, rowname is gene
    # REF: colname is cell type, rowname is gene 
    # COMBAT: use ComBat to do batch-effect correction (default is TRUE)
    # WEIGHT: use gene's variance as weight (default is TRUE)
    
    # mydelia$out: colname is sample name, rowname is cell type
    
    



