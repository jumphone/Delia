
.norm_exp<-function(x){
    y=x
    y[which(x<0)]=0
    y=x/sum(x)
    y=y*1000000
    return(y)
    }

.norm_one <- function(x){
    y=scale(x) 
    y=pnorm(y)
    y=y/sum(y)

    return(y)
    }

.simple_combine <- function(exp_sc_mat1, exp_sc_mat2){    
    exp_sc_mat=exp_sc_mat1
    exp_ref_mat=exp_sc_mat2
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    OUT=list()
    OUT$exp_sc_mat1=exp_sc_mat
    OUT$exp_sc_mat2=exp_ref_mat
    OUT$combine=cbind(exp_sc_mat,exp_ref_mat)
    return(OUT)
    }


.generate_ref <- function(exp_sc_mat, TAG, min_cell=1, refnames=FALSE){
    NewRef=c()
    TAG=as.matrix(TAG)
    if(ncol(TAG)==1){TAG=cbind(TAG,TAG)}
    
    TAG[,2]=as.character(TAG[,2])
    
    if(refnames==FALSE){
        refnames=names(table(TAG[,2]))}
        else{refnames=refnames}
    outnames=c()
    for(one in refnames){
        this_col=which(TAG[,2]==one)
        if(length(this_col)>= min_cell){
            outnames=c(outnames,one)
            if(length(this_col) >1){
                #this_new_ref=apply(exp_sc_mat[,this_col],1,sum)
                this_new_ref=apply(exp_sc_mat[,this_col],1,mean)
                }
                else{this_new_ref = exp_sc_mat[,this_col]}
            NewRef=cbind(NewRef,this_new_ref)
            }
        }
    rownames(NewRef)=rownames(exp_sc_mat)
    colnames(NewRef)=outnames
    if(length(NewRef[1,])==1){
        NewRef=cbind(NewRef[,1], NewRef[,1])
        rownames(NewRef)=rownames(exp_sc_mat)
        colnames(NewRef)=c(outnames,outnames)
        }
    NewRef=apply(NewRef, 2, .norm_exp)
    return(NewRef)
    }



.combat <- function(EXP, BATCH){
    library(sva)
    library(limma)
    pheno = data.frame(batch=as.matrix(BATCH))
    edata = EXP
    batch = pheno$batch
    modcombat = model.matrix(~1, data=pheno)
    combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
    return(combat_edata)
    }


Delia <- function(EXP, REF, COMBAT=TRUE, WEIGHT=TRUE){
    ##############################
    print('Start!')
    print(Sys.time())
    ##############################
    ######
    
    REF=REF
    EXP=EXP
    COMBAT=COMBAT
    WEIGHT=WEIGHT
    
    COM=.simple_combine(EXP, REF)$combine
    NCOM=apply(COM,2,.norm_exp)
    rownames(NCOM)=rownames(COM)
    colnames(NCOM)=colnames(COM)
    COM=NCOM

    ##########
    if(COMBAT==TRUE){
        ##############################
        print('Batch correction...')
        ##############################
        BATCH=c(rep('REF',ncol(REF)),rep('EXP',ncol(EXP)))
        COM.combat=.combat(COM, BATCH) 
        COM.orig=COM
        COM=COM.combat
        NCOM=apply(COM,2,.norm_exp)
        rownames(NCOM)=rownames(COM)
        colnames(NCOM)=colnames(COM)
        COM=NCOM
        }

    ############   
    if(WEIGHT==TRUE){
    	##############################
    	print('Weight calculation...')
    	##############################
        SCOM=t(apply(COM,1,scale))
        rownames(SCOM)=rownames(COM)
        colnames(SCOM)=colnames(COM)
        wgt=apply(SCOM, 1, var)
        COM=SCOM
        }
        
    ##################
    ##############################
    print('Deconvelution ... ')
    ##############################
    OUT=c()
    i=1
    while(i<=ncol(EXP)){
        this_com= COM[,c(i,c((ncol(EXP)+1): ncol(COM) ))]
        colnames(this_com)[1]='NOI'
        this_com=as.data.frame(this_com)
        if(WEIGHT==TRUE){
            fit=lm(NOI ~ ., data=this_com, weights=wgt)
            }else{
            fit=lm(NOI ~ ., data=this_com)	
            }
        this_coef=fit$coefficients
        this_ratio=.norm_one(this_coef[c(2:length(this_coef))])
        OUT=cbind(OUT,this_ratio)
        i=i+1
    }
    rownames(OUT)=colnames(REF)
    colnames(OUT)=colnames(EXP)


    ##################################
    RESULT=list()
    RESULT$exp=EXP
    RESULT$ref=REF
    RESULT$combat=COMBAT
    
    RESULT$out=OUT
    
    if(COMBAT==TRUE){
        RESULT$combat.exp=COM.orig
        RESULT$combat.exp=COM.combat
        RESULT$combat.batch=BATCH
        }
    ##############################
    print('Finished!')
    print(Sys.time())
    ##############################
    return(RESULT)
    }


