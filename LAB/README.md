
estimate_ratio_pv=-log(mydelia$pvalue,10)

True.A=true_ratio[1,]
Est.A=estimate_ratio_pv[1,]
PCC=round(cor(True.A,Est.A,method='pearson'),2)

plot(True.A, Est.A, pch=16, main=paste0('PCC=',PCC))

RATIO=estimate_ratio_pv

ALLR=true_ratio
CORMAT=cor(t(RATIO), t(ALLR), method='pearson')
