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
    
    # Result is in "mydelia$out"
    
    



