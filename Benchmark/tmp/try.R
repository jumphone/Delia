




get_comparison_markers <- function(method, data, q, pure_samples, data_type, sig, 
    all_markers = FALSE) {
    
    updt(paste(sig, "Finding Markers"))
    pure <- unlist(pure_samples)
    K <- length(pure_samples)
    n_choose <- floor((1 - q) * ncol(data)/K)
    
    if (method == "Ratio" || method == "Diff" || method == "p.value") {
        marker_list_ratio <- find_markers2(data, pure_samples, data_type = data_type, 
                                           marker_method = tolower(method))
        ## Choose number of markers
        q_ratio <- sapply(1:K, function(i) {
            quantile(marker_list_ratio$V[[i]], q)
        })
        med_q <- median(q_ratio)
        chs <- function(q_ratio, j) {
            qtn <- q_ratio[j]
            n <- max(which(marker_list_ratio$V[[j]] >= qtn), 1)
            return(n)
        }
        n_choose_ratio <- sapply(1:K, function(a) {
            chs(q_ratio, a)
        })
        marks_ratio <- marker_list_ratio$L
        mrks <- lapply(1:K, function(i) {
            marks_ratio[[i]][1:n_choose_ratio[i]]
        })
        if (all_markers) {
            mrks <- marks_ratio
        }
    } else if (method == "p-value") {
        factors <- unlist(lapply(1:K, function(x) rep(names(pure_samples)[x], length(pure_samples[[x]]))))
        ufactors <- names(pure_samples)
        markers_ged <- extractMarkers(t(data)[, pure], factors, method = "Abbas", 
            log = FALSE)
        for (nm in ufactors) {
            markers_ged[[nm]] <- markers_ged[[nm]][is.finite(markers_ged[[nm]])]
        }
        
        if (all(lengths(markers_ged) == 0)) {
            marker_mtx <- data.frame(markerScoreAbbas(t(data)[, pure], data = factor(factors), 
                log = FALSE))
            markers_ged <- split(marker_mtx$dm2, marker_mtx$top)
            L <- split(rownames(marker_mtx), marker_mtx$top)
            for (i in 1:length(markers_ged)) {
                names(markers_ged[[i]]) <- L[[i]]
            }
            markers_ged <- lapply(markers_ged, sort, decreasing = TRUE)
            names(markers_ged) <- sapply(factors, toString)
        }
        qs <- lapply(ufactors, function(x) {
            quantile(markers_ged[[x]], 1 - q)
        })
        mrks <- lapply(1:K, function(i) {
            which(colnames(data) %in% names(which(markers_ged[[ufactors[i]]] <= qs[i])))
        })
    }
    mrks <- lapply(mrks, function(x) x[1:min(n_choose, length(x))])
    names(mrks) <- names(pure_samples)
    return(mrks)
}





setwd('F:/Delia/COMPARE/dtangle.data_0.3.6.tar/dtangle.data/data')
source('https://raw.githubusercontent.com/jumphone/Delia/master/Delia.R')

load('abbas.rda')
DATA=abbas


EXP=t(DATA$data$log)

REF.FLAG=apply(DATA$annotation$pure,1,sum)
which(REF.FLAG==1)

REF.TAG=rep(NA,length(REF.FLAG))
i=1
while(i<=nrow(DATA$annotation$pure)){
    if(REF.FLAG[i]==1){
        this_tag=colnames(DATA$annotation$pure)[which(DATA$annotation$pure[i,]==TRUE)]
        REF.TAG[i]=this_tag
        }

    i=i+1}

REF.INDEX=which(REF.FLAG==1)

REF=.generate_ref(EXP[,REF.INDEX], REF.TAG[REF.INDEX])


MIX.INDEX=which(REF.FLAG==0)
MIX=EXP[,MIX.INDEX]


mydelia <- Delia(MIX, REF, COMBAT=TRUE, METHOD='glm', SHOW=FALSE, RANK=FALSE)     

PP.TRUE=DATA$annotation$mixture[MIX.INDEX,]
PP.DELIA=t(mydelia$out)

PP.TRUE=PP.TRUE[,order(colnames(PP.TRUE))]
PP.DELIA=PP.DELIA[,order(colnames(PP.DELIA))]











