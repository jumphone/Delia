source('Delia.R')
sc_exp_mat <- .readTable(PATH='SingleCellREF/sc_exp_mat.txt', SEP='\t')
sc_exp_mat=log(sc_exp_mat+1,10)
sc_exp_mat=apply(sc_exp_mat, 2, .norm_exp)

tag <- read.table('SingleCellREF/tag.txt', sep='\t', header=TRUE,row.names=1)
tag <- as.character(tag[,1])

EXP <- readRDS('EXP.RDS')  # This matrix has been log-normalized already. 
true_ratio <- readRDS('TRUE_RATIO.RDS')
REF <- .generate_ref(sc_exp_mat, tag)



mydelia <- Delia(EXP, REF, COMBAT=TRUE, METHOD='opt', SHOW=FALSE)      










N=1000

PseudoData<-data.frame(x1=runif(N, -5, 5),
                        x2=runif(N, -5, 5), 
                        x3=runif(N, -5, 5)) 
 

PseudoData<-within(PseudoData, { 
   y1=1+x1+0.3*x2+0.7*x3+rnorm(N, 0, 3) 
   y2=1+x1+0.5*x2-0.5*x3+rnorm(N, 0, 3) 
 })


Data1<-with(PseudoData,data.frame(y=y1, x1=x1, x2=x2, x3=x3)) 
Data2<-with(PseudoData,data.frame(y=y2, x1=x1, x2=x2, x3=x3)) 
 


eval.e="e<-expression((y-(c1+x1+c2*x2+c3*x3))^2)"
eval.foo="foo<-deriv(e, nam=c('c1','c2','c3'))"






DATA=Data1




lm.optim <- function(DATA){

    CN=colnames(DATA)
    VAR_NAME=CN[2:length(CN)]
    COEF_NAME=paste0('C', c(1:(length(CN)-1)))
    ALL_COEF=c('C0',COEF_NAME)
    ALL_COEF_E=paste(ALL_COEF,collapse ="','")
    E=paste(paste0(COEF_NAME,'*',VAR_NAME),collapse="+")

    eval.e=paste0("e<-expression((",CN[1],"-(C0+",E,"))^2)")
    eval.foo=paste0("foo<-deriv(e, nam=c('", ALL_COEF_E, "'))")

    eval(parse(text = eval.e))
    eval(parse(text = eval.foo))

    LW=c(-Inf, rep(0,length(ALL_COEF)-1))
    UP=c(Inf, rep(1,length(ALL_COEF)-1))
    PAR=rep(0.5, length(ALL_COEF))
    names(PAR)=ALL_COEF

    objfun<-function(coefs, data) { 
       return(sum(eval(foo,env=c(as.list(coefs), as.list(data))))) 
     } 
 
    objgrad<-function(coefs, data) { 
        return(apply(attr(eval(foo,env=c(as.list(coefs), as.list(data))), 
                        "gradient"),2,sum)) 
     } 

    DATA=as.data.frame(DATA)

    D1.bound<-optim(par=PAR, 
                fn=objfun, 
                gr=objgrad, 
                data=DATA, 
                method="L-BFGS-B", 
                lower=LW, 
                upper=UP)


    OUT=D1.bound$par
    names(OUT)[2:length(OUT)]=VAR_NAME
    return(OUT)
    }



OUT=lm.optim(DATA)











.beta_opt <- function(mydelia, N=30, s2=3){
    N=N
    s2=s2
    s1.lst=c(1:N) *s2 /N    
    COEF=mydelia$coef
    COM=mydelia$com
    NREF=COM[,(1+ncol(mydelia$exp)):ncol(COM)]
    NEXP=COM[,c(1:ncol(mydelia$exp))]
    SUM=c()
    I=c()

    i=1
    while(i<=N){
        this_s1=s1.lst[i]
        this_s2=s2   
        OUT=apply(COEF,2,.beta_one, s1=this_s1,s2=this_s2)
        predict.EXP=NREF %*% OUT
        this_index=1
        this_cor.lst=c()
        while(this_index<=ncol(NEXP)){
            this_cor=cor(NEXP[,this_index],predict.EXP[,this_index],method='spearman')
            this_cor.lst=c(this_cor.lst,this_cor)
            this_index=this_index+1}
        this_sum=sum(this_cor.lst)
        SUM=c(SUM, this_sum)
        I=c(I,i)
        print(i)
        i=i+1
        }

    s1.opt=s1.lst[I[which(SUM==max(SUM))]]
    s2.opt=s2
    OUT=apply(COEF,2,.beta_one, s1=s1.opt,s2=s2.opt)
    rownames(OUT)=rownames(COEF)
    colnames(OUT)=colnames(COEF)
    RESULT=list()
    RESULT$out=OUT
    RESULT$s1=s1.opt
    RESULT$s2=s2.opt
    return(RESULT)
    }


OUT=.beta_opt(mydelia , N=30,s2=3)
OOO=OUT$out

cor(OOO[1,],true_ratio[1,])

    #print(s1.opt)
    #print(s2.opt)\
hist(rbeta(10000,shape1=OUT$s1,shape2=OUT$s2))

OUT=apply(COEF,2,.beta_one, s1=1,s2=1)
cor(OUT[1,],true_ratio[1,])


OUT=apply(COEF,2,.norm_one)
cor(OUT[1,],true_ratio[1,])

















setwd('F:/Delia')
# Load Data
REF=readRDS('./RDS/REF.RDS')
ALLR=readRDS('./RDS/ALLR.RDS')
EXP=readRDS('./RDS/EXP.RDS')
REXP=readRDS('./RDS/REXP.RDS')
SC.REF=readRDS('./RDS/SC.REF.RDS')
V.SC.REF=readRDS('./RDS/V.SC.REF.RDS')
V.SC.REF_scmat=readRDS('./RDS/V.SC.REF_scmat.RDS')
################
VAR=apply(SC.REF, 1, var)
VG=which(rank(-VAR) <= 2000  )
########################

library(deconvSeq)
############################


sc_mat=read.table('Zeisel_exp_sc_mat.txt',sep='\t',header=T,row.names=1)
sc_mat[which(is.na(sc_mat))]
SCMAT=sc_mat[VG,]

colnames(SCMAT)=paste0(colnames(SCMAT),'_',c(1:ncol(SCMAT)))
tissue.sc = colnames(V.SC.REF_scmat)
design.sc = model.matrix(~-1+as.factor(tissue.sc))
colnames(design.sc) = levels(as.factor(tissue.sc))
rownames(design.sc) = names(tissue.sc)
dge.sc = getdge(SCMAT, design.sc, ncpm.min=1, nsamp.min=4, method="bin.loess")
b0.sc = getb0.rnaseq(dge.sc, design.sc, ncpm.min=1, nsamp.min=4)
##########################################


source('Delia.R')
COM=.simple_combine(REXP, SCMAT)$com
BATCH=c(rep('EXP',ncol(REXP)), rep('REF',ncol(SCMAT)))
COM.combat=.combat(COM, BATCH) 


####################################################
dge.tissue = getdge(REXP,design=NULL, ncpm.min=1,nsamp.min=4)
####################################################
resultx1_s.sc = getx1.rnaseq(NB0=50, b0.sc, dge.tissue)






# Load Package 
library(sva)
library(limma)
################

source('Delia.R')
# Delia with all genes, 5 seconds
mydelia = Delia(REXP, SC.REF) 
RATIO=mydelia$out
#######
saveRDS(RATIO, file='./RESULT/Delia_all.RDS')
#######
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
pdf('./RESULT/Delia_all.pdf',width=6, height=6)
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
dev.off()
########################



RATIO=readRDS(file='./RESULT/Delia_all.RDS')
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]

plot(RATIO[4,],ALLR[6,],pch=16)


# Delia with variable genes, 1 second
mydelia = Delia(REXP, V.SC.REF) 
RATIO=mydelia$out
#######
saveRDS(RATIO, file='./RESULT/Delia_var.RDS')
#######
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
pdf('./RESULT/Delia_var.pdf',width=6, height=6)
library('gplots')
heatmap.2(CORMAT,scale=c("none"),dendrogram='none',Rowv=F,Colv=F,cellnote=round(CORMAT,2),notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),margins=c(10,10))
dev.off()
#########################

RATIO=readRDS(file='./RESULT/Delia_var.RDS')
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
CORMAT[1,1]+CORMAT[2,2]+CORMAT[3,4]+CORMAT[3,5]+CORMAT[4,6]+CORMAT[5,7]








