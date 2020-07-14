### PCA分析，观察组间差异是否明显

```R
load("gdc.Rdata")
library(FactoMineR)
library(factoextra)
dat <- log(expr+1)
dat.pca <- PCA(t(dat),graph = F)
pca.plot <- fviz_pca_ind(dat.pca,geom.ind = 'point',col.ind = group_list,
                         addEllipses = T,legend.title = 'Groups')
```

