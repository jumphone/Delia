
print("Welcome to use Delia(v0.1.2) ! ")

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



.beta_one <- function(x, s1=0.1,s2=10){
    s1=s1
    s2=s2
    y=scale(x) 
    y=pnorm(y)
    y=qbeta(y, shape1=s1,shape2=s2)
    y=y/sum(y)

    return(y)
    }



lm.optim <- function(DATA){

    CN=colnames(DATA)
    VAR_NAME=CN[2:length(CN)]
    COEF_NAME=paste0('C', c(1:(length(CN)-1)))
    ALL_COEF=c('C0',COEF_NAME)
    ALL_COEF_E=paste(ALL_COEF,collapse ="','")
    E=paste(paste0(COEF_NAME,'*',VAR_NAME),collapse="+")
    ##################
    eval.e=paste0("e<-expression((",CN[1],"-(C0+",E,"))^2)")
    eval.foo=paste0("foo<-deriv(e, nam=c('", ALL_COEF_E, "'))")
    ##################
    eval(parse(text = eval.e))
    eval(parse(text = eval.foo))
    ##################
    LW=c(0, rep(0,length(ALL_COEF)-1))
    UP=c(1, rep(1,length(ALL_COEF)-1))
    PAR=rep(0.5, length(ALL_COEF))
    names(PAR)=ALL_COEF
    ###################
    objfun<-function(coefs, data) { 
       return(sum(eval(foo,env=c(as.list(coefs), as.list(data))))) 
     } 
    ###################
    objgrad<-function(coefs, data) { 
        return(apply(attr(eval(foo,env=c(as.list(coefs), as.list(data))), 
                        "gradient"),2,sum)) 
     } 
    ###################
    DATA=as.data.frame(DATA)
    D1.bound<-optim(par=PAR, 
                fn=objfun, 
                gr=objgrad, 
                data=DATA, 
                method="L-BFGS-B", 
                lower=LW, 
                upper=UP)
    ####################
    OUT=D1.bound$par
    names(OUT)[2:length(OUT)]=VAR_NAME
    return(OUT)
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




Delia <- function(EXP, REF, COMBAT=TRUE, RANK=FALSE, SHOW=FALSE, METHOD='lm', PCV=0.95){
    ##############################
    print('Start!')
    print(Sys.time())
    ##############################   
    REF=REF
    EXP=EXP
    #############################
    header.REF=colnames(REF)
    header.EXP=colnames(EXP)
    colnames(REF)=paste0('REF.',c(1:ncol(REF)))
    colnames(EXP)=paste0('EXP.',c(1:ncol(EXP)))
    ############################
    COMBAT=COMBAT
    METHOD=METHOD
    RANK=RANK
    ##############################
    # Check package
    if(METHOD=='pcr'){
        if(!'pls' %in% installed.packages()[,1]){
            print("Please install 'pls' for the principal components regression.")
            print("install.packages('pls')")
            return(NULL)}
        library(pls)
        }
    if(METHOD=='rlm'){
        if(!'MASS' %in% installed.packages()[,1]){
            print("Please install 'MASS' for the robust linear regression.")
            print("install.packages('MASS')")
            return(NULL)}
        library(MASS)
        }
    
    ###############################
    if(SHOW==TRUE){
        if(!'tcltk2' %in% installed.packages()[,1]){
            print("Please install 'tcltk2' for the progress bar.")
            print("install.packages('tcltk2')")
            return(NULL)}
        library(tcltk2)}
    ###############################
    PCV=PCV
    ##############################
    if(METHOD=='opt'){
        print("Because METHOD is 'opt', RANK is changed to TRUE.")
        RANK=TRUE}
    ##############################
    #SCALE=TRUE
    #RANK=TRUE
    
    if(RANK==TRUE){
        SCALE=FALSE
    }else{
        SCALE=TRUE
    }
    
    ###############################
    COM=.simple_combine(EXP, REF)$combine
    NCOM=apply(COM,2,.norm_exp)
    rownames(NCOM)=rownames(COM)
    colnames(NCOM)=colnames(COM)
    ###########################
    VAR=apply(NCOM, 1, var)
    NCOM=NCOM[which(VAR>0),]
    ##########################
    COM=NCOM
    ##########
    if(COMBAT==TRUE){
        ##############################
        print('Batch correction...')
        ##############################
        BATCH=c(rep('REF',ncol(REF)),rep('EXP',ncol(EXP)))
        COM.combat=.combat(COM, BATCH)         
        COM=COM.combat
        NCOM=apply(COM,2,.norm_exp)
        rownames(NCOM)=rownames(COM)
        colnames(NCOM)=colnames(COM)
        COM=NCOM
        }
    ############
    #####################
    COM.orig=COM
    ##################
    if(RANK==TRUE){
        ##############################
        print('Rank-normalization...')
        ##############################    
        COM=apply(COM.orig,2,rank)
        }
    ##################
    ############  
    if(SCALE==TRUE){
        ##############################
        print('Standardization...')
        ##############################    
        SCOM=t(apply(COM,1,scale))
        rownames(SCOM)=rownames(COM)
        colnames(SCOM)=colnames(COM)
        COM=SCOM}
    ##################
    ##############################
    print('Deconvolution ... ')
    ##############################
    OUT=c()
    C=c()
    PCN=c()
    OPT=c()
    NN=ncol(EXP)
    ###################
    if(SHOW==TRUE){
        pb = tkProgressBar('Delia',"/", 0, NN)
        }
    #####################
    i=1
    while(i<=ncol(EXP)){
        this_com= COM[,c(i,c((ncol(EXP)+1): ncol(COM) ))]
        colnames(this_com)[1]='NOI'
        this_com=as.data.frame(this_com)
        ##############################
        if(METHOD=='pcr'){
            fit=pcr(NOI~., data = this_com, scale = TRUE, validation = "CV")
            Xvar=fit$Xvar/fit$Xtotvar
            used_pc=1
            while(sum(Xvar[1:used_pc])<PCV & used_pc<ncol(REF) ){used_pc=used_pc+1}
            PCN=c(PCN,used_pc)
            ####################
            this_coef=fit$coefficients[,,used_pc]
            this_ratio=.norm_one(this_coef)
            }
        ##############################
        if(METHOD=='lm'){
            fit=lm(NOI ~ ., data=this_com)  
            ####################
            this_coef=fit$coefficients[c(2:(ncol(REF)+1))]
            this_ratio=.norm_one(this_coef)
            }
        ##############################
        if(METHOD=='rlm'){
            fit=rlm(NOI ~ ., data=this_com) 
            #################### 
            this_coef=fit$coefficients[c(2:(ncol(REF)+1))]
            this_ratio=.norm_one(this_coef)
            }
        ############################
        ############################
        if(METHOD=='opt'){       
            opt.out=lm.optim(this_com)
            OPT=cbind(OPT,opt.out)
            ################
            this_coef=opt.out[c(2:(length(this_com)))]
            this_ratio=(this_coef/sum(this_coef))
            }
        ############################
        ############################
        
        #############################        
        OUT=cbind(OUT,this_ratio)
        C=cbind(C, this_coef)
        ############################
        if(SHOW==TRUE){
            info = sprintf(paste0("%d / ",NN), round(i))
            setTkProgressBar(pb, i, sprintf('Delia (%s)',info), info)
            }      
        ###################
        i=i+1
    }
    rownames(OUT)=header.REF
    colnames(OUT)=header.EXP
    rownames(C)=header.REF
    colnames(C)=header.EXP
    
    ##################################
    RESULT=list()
    RESULT$combat=COMBAT
    RESULT$out=OUT
    RESULT$coef=C
    RESULT$pcn=PCN  
    RESULT$pcv=PCV
    RESULT$method=METHOD
    RESULT$rank=RANK
    ######################
    if(SHOW==TRUE){
        close(pb)
        }        
    ###############################
    print('Finished!')
    print(Sys.time())
    ##############################
    return(RESULT)
    }





.readTable <- function(PATH, SEP='\t', UP=FALSE){
    SEP=SEP
    PATH=PATH
    UP=UP
    DATA=read.table(file=PATH,sep=SEP,header=TRUE,row.names=NULL)
    DATA=apply(DATA,2,as.character)
    ###########
    if(UP==TRUE){DATA[,1]=toupper(DATA[,1])}
    ###########
    TAB=table(as.character(DATA[,1]))
    UNIQ=names(TAB)[which(TAB==1)]
    DATA=DATA[which(DATA[,1] %in% UNIQ),]
    RN=DATA[,1]
    DATA=DATA[,c(2:ncol(DATA))]
    DATA=apply(DATA,2,as.numeric)
    rownames(DATA)=RN
    return(DATA)
    }


.writeTable <- function(DATA, PATH, SEP='\t',TL='OUTPUT'){
    DATA=DATA
    PATH=PATH
    SEP=SEP
    TL=TL
    OUT=cbind(rownames(DATA),DATA)
    colnames(OUT)[1]=TL
    write.table(OUT, file=PATH,sep=SEP,row.names=FALSE,col.names=TRUE,quote=FALSE)
    }
    
    
.scaleCoef <- function(COEF){
    COEF=COEF
    SCOEF=apply(COEF,2,scale)
    SCOEF=t(apply(SCOEF,1,scale))
    rownames(SCOEF)=rownames(COEF)
    colnames(SCOEF)=colnames(COEF)
    return(SCOEF)
    }


.scaleGlobal <- function(DATA){
    DATA=DATA#mydelia$coef
    OUT=matrix(scale(as.numeric(DATA)), ncol=ncol(DATA),nrow=nrow(DATA))
    #OUT=pnorm(OUT)
    rownames(OUT)=rownames(DATA)
    colnames(OUT)=colnames(DATA)
    return(OUT)
    }


