Delia <- function(EXP, REF, COMBAT=TRUE, SHOW=FALSE, METHOD='lm', PCV=0.95){
    ##############################
    print('Start!')
    print(Sys.time())
    ##############################   
    REF=REF
    EXP=EXP
    COMBAT=COMBAT
    METHOD=METHOD
    if(METHOD=='pcr'){
        if(!'pls' %in% installed.packages()[,1]){
            print("Please install 'pls' for the PCR.")
            print("install.packages('pls')")
            return(NULL)}
        library(pls)
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
    SCALE=TRUE
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
    COM.orig=COM
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
    NN=ncol(EXP)
    ###################
    if(SHOW==TRUE){
        pb = tkProgressBar('Progress',"Finished %", 0, 100)
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
            this_coef=fit$coefficients[,,used_pc]
            }
        ##############################
        if(METHOD=='lm'){
            fit=lm(NOI ~ ., data=this_com)  
            this_coef=fit$coefficients[c(2:(ncol(REF)+1))]
            }
        ##############################
        if(METHOD=='rlm'){
            library(MASS)
            fit=rlm(NOI ~ ., data=this_com)  
            this_coef=fit$coefficients[c(2:(ncol(REF)+1))]
            }

        ############################
        this_ratio=.norm_one(this_coef)
        #############################        
        OUT=cbind(OUT,this_ratio)
        C=cbind(C, this_coef)
        ############################
        if(SHOW==TRUE){
            info = sprintf("Finished %d%%", round(i*100/NN))
            setTkProgressBar(pb, i*100/NN, sprintf('Progress (%s)',info), info)
            }      
        ###################
        i=i+1
    }
    rownames(OUT)=colnames(REF)
    colnames(OUT)=colnames(EXP)
    rownames(C)=colnames(REF)
    colnames(C)=colnames(EXP)
    
    ##################################
    RESULT=list()
    RESULT$exp=EXP
    RESULT$ref=REF
    RESULT$com=COM.orig
    RESULT$combat=COMBAT
    ######################  
    RESULT$out=OUT
    RESULT$coef=C
    RESULT$pcn=PCN
    ######################
    if(COMBAT==TRUE){
        RESULT$combat.exp=COM.combat
        RESULT$combat.batch=BATCH
        }
    ##############################
    if(SHOW==TRUE){
        close(pb)
        }        
    ###############################
    print('Finished!')
    print(Sys.time())
    ##############################
    return(RESULT)
    }
