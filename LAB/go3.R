


EXP=.readTable('BT319_714_mix.txt',SEP='\t')
REF=.readTable('CMAT.txt.rmCC.txt',SEP='\t')

REF=REF[,which(colnames(REF)!='X30')]

source('Delia.R')
mydelia <- Delia(EXP, REF, COMBAT=TRUE, METHOD='opt',SHOW=TRUE)  



#CCC=.simple_combine(EXP,REF)
#CCC$exp_sc_mat1
#CCC$exp_sc_mat2


cor(CCC$exp_sc_mat1,CCC$exp_sc_mat2)

show_ratio <-  mydelia$coef

library('gplots')
 
heatmap.2(t(show_ratio),scale=c("none"), dendrogram='none',
    Rowv=F,Colv=F,cellnote=round(t(show_ratio),2), notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),
    margins=c(10,10))





EXP=.readTable('BT309_G3.txt',SEP='\t')
REF=.readTable('CMAT.txt.rmCC.txt',SEP='\t')
REF=REF[,which(colnames(REF)!='X30')]


source('Delia.R')
mydelia <- Delia(EXP, REF, COMBAT=TRUE, METHOD='opt',SHOW=TRUE)  

show_ratio <-  mydelia$coef

#show_ratio <-  .scaleGlobal(mydelia$coef)

library('gplots')
heatmap.2(t(show_ratio),scale=c("none"), dendrogram='none',
    Rowv=F,Colv=F,cellnote=round(t(show_ratio),2), notecol='black',
    trace='none',col=colorRampPalette(c('blue','royalblue','grey80','indianred','red')),
    margins=c(10,10))


EXP=.readTable('BT308_G4.txt',SEP='\t')
REF=.readTable('CMAT.txt.rmCC.txt',SEP='\t')
REF=REF[,which(colnames(REF)!='X30')]


mydelia <- Delia(EXP, REF, COMBAT=TRUE, METHOD='opt',SHOW=TRUE)  


show_ratio <-  mydelia$coef

library('gplots')
 
heatmap.2(t(show_ratio),scale=c("none"), dendrogram='none',
    Rowv=F,Colv=F,cellnote=round(t(show_ratio),2), notecol='black',
    trace='none',col=colorRampPalette(c('royalblue','grey80','indianred')),
    margins=c(10,10))

