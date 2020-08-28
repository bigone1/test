an R/Bioconductor package for downloading, organizing, and integrative analyzing lncRNA, mRNA, and miRNA data in GDC    https://github.com/Jialab-UCR/GDCRNATools

```R
library(GDCRNATools)
#mRNA的表达矩阵
data(rnaCounts)
#临床信息
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'RNAseq',
                                   write.meta = FALSE)
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)
#过滤
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)
#差异分析
DEGAll <- gdcDEAnalysis(counts     = rnaCounts, 
                        group      = metaMatrix.RNA$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma')
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')
#两个基因的相关性图
gdcCorPlot(gene1    = 'ENSG00000003989', 
           gene2    = 'ENSG00000004799', 
           rna.expr = rnaExpr, 
           metadata = metaMatrix.RNA)
#生存分析：两种方法---coxph、KM
survOutput1 <- gdcSurvivalAnalysis(gene     = rownames(dePC), 
                                  method   = 'coxph', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA)
survOutput2 <- gdcSurvivalAnalysis(gene     = rownames(dePC), 
                                  method   = 'KM', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA, 
                                  sep      = 'median')
gdcKMPlot(gene     = 'ENSG00000003989',
          rna.expr = rnaExpr,
          metadata = metaMatrix.RNA,
          sep      = 'median')
```

