```R
library(estimate)
OvarianCancerExpr <- system.file('extdata','sample_input.txt',package = 'estimate')
filterCommonGenes(input.f = OvarianCancerExpr,output.f = 'genes.gct',id='GeneSymbol')
estimateScore(input.ds = 'genes.gct',output.ds = 'score.gct')#参数platform，默认为'affymetrix',还有'illumina','agilent'
scores=read.table("score.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
write.csv(scores,file = 'scores.csv')
```
