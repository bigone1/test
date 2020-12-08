## READ

```R
read <- data.table::fread('lianxi/PMID32798715/READ/TCGA-READ.htseq_fpkm.tsv')
surv <- data.table::fread('lianxi/PMID32798715/READ/TCGA-READ.survival.tsv')
surv <- surv[,-3]
pheno <- data.table::fread('lianxi/PMID32798715/READ/TCGA-READ.GDC_phenotype.tsv')
pheno <- pheno[,c("submitter_id.samples","age_at_initial_pathologic_diagnosis","gender.demographic","pathologic_M",
                  "pathologic_N","pathologic_T","tumor_stage.diagnoses")]
pheno <- na.omit(pheno)
pheno <- pheno[which(pheno$pathologic_T != "TX"),]
pheno[which(pheno$pathologic_T == "T1a"),6] <- "T1"
pheno[which(pheno$pathologic_T == "T1b"),6] <- "T1"
pheno[which(pheno$pathologic_T == "T2a"),6] <- "T2"
pheno[which(pheno$pathologic_T == "T2b"),6] <- "T2"
pheno[which(pheno$pathologic_M == "M1a"),4] <- "M1"
pheno[which(pheno$pathologic_M == "M1b"),4] <- "M1"
pheno[which(pheno$pathologic_N == "N1b"),5] <- "N1"
pheno[which(pheno$pathologic_N == "N1a"),5] <- "N1"
pheno[which(pheno$pathologic_N == "N1c"),5] <- "N1"
pheno[which(pheno$pathologic_N == "N2a"),5] <- "N2"
pheno[which(pheno$pathologic_N == "N2b"),5] <- "N2"
pheno[which(pheno$pathologic_T == "T4a"),6] <- "T4"
pheno[which(pheno$pathologic_T == "T4b"),6] <- "T4"
pheno <- pheno[which(pheno$pathologic_N != "NX"),]
pheno <- pheno[which(pheno$pathologic_M != ""),]
pheno <- pheno[which(pheno$tumor_stage.diagnoses != "not reported"),]

pheno[which(pheno$tumor_stage.diagnoses == "stage i"),7] <- "I"
pheno[which(pheno$tumor_stage.diagnoses == "stage ia"),7] <- "I"
pheno[which(pheno$tumor_stage.diagnoses == "stage ib"),7] <- "I"

pheno[which(pheno$tumor_stage.diagnoses == "stage ii"),7] <- "II"
pheno[which(pheno$tumor_stage.diagnoses == "stage iia"),7] <- "II"
pheno[which(pheno$tumor_stage.diagnoses == "stage iib"),7] <- "II"
pheno[which(pheno$tumor_stage.diagnoses == "stage iic"),7] <- "II"
pheno[which(pheno$tumor_stage.diagnoses == "stage iii"),7] <- "III"
pheno[which(pheno$tumor_stage.diagnoses == "stage iiia"),7] <- "III"
pheno[which(pheno$tumor_stage.diagnoses == "stage iiib"),7] <- "III"
pheno[which(pheno$tumor_stage.diagnoses == "stage iiic"),7] <- "III"
pheno[which(pheno$tumor_stage.diagnoses == "stage iv"),7] <- "IV"
pheno[which(pheno$tumor_stage.diagnoses == "stage iva"),7] <- "IV"
exp_sample <- colnames(read)[-1]
surv_sample <- surv$sample
pheno_sample <- pheno$submitter_id.samples
jiaoji <- intersect(intersect(exp_sample,surv_sample),pheno_sample)

pheno1 <- pheno[which(pheno$submitter_id.samples %in% jiaoji),]
colnames(pheno1) <- c("sample_id","age","gender","M","N","T","stage")
read <- tibble::column_to_rownames(read,"Ensembl_ID")
read1 <- read[,jiaoji] 

surv1 <- surv[which(surv$sample %in% jiaoji),]

irg <- data.table::fread('D:/linshiwenjian/r-daima/免疫基因.txt')
irg <- irg[!duplicated(irg$Symbol),]
idmapping <- read.table(file = "lianxi/PMID32798715/gencode.v22.annotation.gene.probeMap", 
                        sep = "\t", header = T, stringsAsFactors = F)
geneid <- data.frame(id = rownames(read1), stringsAsFactors = F)
geneid2symbol <- dplyr::left_join(geneid, idmapping, by = "id")
read1$symbol <- geneid2symbol$gene
read1 <- read1[which(read1$symbol %in% irg$Symbol),]
read1 <- read1[order(read1$symbol),]
read2 <- aggregate(read1[,c(1:(ncol(read1)-1))], 
                   by = list(read1$symbol), mean)
rownames(read2) <- read2[,1]
read2 <- read2[,-1]
irgps <- read.table("lianxi/PMID32798715/READ/IRGPs.txt")
lieming <- colnames(irgps)
lieming1 <- unlist(lapply(lieming,gsub,pattern="[.]",replacement="-"))
colnames(irgps) <- lieming1
irgps1 <- irgps[,colnames(read2)]
a <- rownames(irgps1)
for(i in a){
  irgps1[i,] <- (read2[unlist(strsplit(i,split='[|]'))[1],]+read2[unlist(strsplit(i,split='[|]'))[2],])/2
}
save(read2,file = "lianxi/PMID32798715/READ/read_exp.RData")

```

## STAD

```R
stad <- data.table::fread('lianxi/PMID32798715/STAD/TCGA-STAD.htseq_fpkm.tsv')
surv <- data.table::fread('lianxi/PMID32798715/STAD/TCGA-STAD.survival.tsv')
surv <- surv[,-3]
pheno <- data.table::fread('lianxi/PMID32798715/STAD/TCGA-STAD.GDC_phenotype.tsv')
pheno <- pheno[,c("submitter_id.samples","age_at_initial_pathologic_diagnosis","gender.demographic","pathologic_M",
                  "pathologic_N","pathologic_T","tumor_stage.diagnoses")]
#pheno <- na.omit(pheno)
pheno <- pheno[which(pheno$pathologic_T != "TX"),]
pheno <- pheno[which(pheno$pathologic_N != "NX"),]
table(pheno$pathologic_M)
pheno <- pheno[which(pheno$pathologic_M != ""),]
table(pheno$pathologic_N)
pheno <- pheno[which(pheno$pathologic_N != ""),]
table(pheno$pathologic_T)
table(pheno$tumor_stage.diagnoses)
pheno <- pheno[which(pheno$tumor_stage.diagnoses != "not reported"),]
pheno[which(pheno$pathologic_N == "N3a"),5] <- "N3"
pheno[which(pheno$pathologic_N == "N3b"),5] <- "N3"
pheno[which(pheno$pathologic_T == "T1a"),6] <- "T1"
pheno[which(pheno$pathologic_T == "T1b"),6] <- "T1"
pheno[which(pheno$pathologic_T == "T2a"),6] <- "T2"
pheno[which(pheno$pathologic_T == "T2b"),6] <- "T2"
pheno[which(pheno$pathologic_M == "M1a"),4] <- "M1"
pheno[which(pheno$pathologic_M == "M1b"),4] <- "M1"
pheno[which(pheno$pathologic_N == "N1b"),5] <- "N1"
pheno[which(pheno$pathologic_N == "N1a"),5] <- "N1"
pheno[which(pheno$pathologic_N == "N1c"),5] <- "N1"
pheno[which(pheno$pathologic_N == "N2a"),5] <- "N2"
pheno[which(pheno$pathologic_N == "N2b"),5] <- "N2"
pheno[which(pheno$pathologic_T == "T4a"),6] <- "T4"
pheno[which(pheno$pathologic_T == "T4b"),6] <- "T4"

pheno[which(pheno$tumor_stage.diagnoses == "stage i"),7] <- "I"
pheno[which(pheno$tumor_stage.diagnoses == "stage ia"),7] <- "I"
pheno[which(pheno$tumor_stage.diagnoses == "stage ib"),7] <- "I"

pheno[which(pheno$tumor_stage.diagnoses == "stage ii"),7] <- "II"
pheno[which(pheno$tumor_stage.diagnoses == "stage iia"),7] <- "II"
pheno[which(pheno$tumor_stage.diagnoses == "stage iib"),7] <- "II"
pheno[which(pheno$tumor_stage.diagnoses == "stage iic"),7] <- "II"
pheno[which(pheno$tumor_stage.diagnoses == "stage iii"),7] <- "III"
pheno[which(pheno$tumor_stage.diagnoses == "stage iiia"),7] <- "III"
pheno[which(pheno$tumor_stage.diagnoses == "stage iiib"),7] <- "III"
pheno[which(pheno$tumor_stage.diagnoses == "stage iiic"),7] <- "III"
pheno[which(pheno$tumor_stage.diagnoses == "stage iv"),7] <- "IV"
pheno[which(pheno$tumor_stage.diagnoses == "stage iva"),7] <- "IV"
exp_sample <- colnames(stad)[-1]
surv_sample <- surv$sample
pheno_sample <- pheno$submitter_id.samples
jiaoji <- intersect(intersect(exp_sample,surv_sample),pheno_sample)

pheno1 <- pheno[which(pheno$submitter_id.samples %in% jiaoji),]
colnames(pheno1) <- c("sample_id","age","gender","M","N","T","stage")
stad <- tibble::column_to_rownames(stad,"Ensembl_ID")
stad1 <- stad[,jiaoji] 

surv1 <- surv[which(surv$sample %in% jiaoji),]

irg <- data.table::fread('D:/linshiwenjian/r-daima/免疫基因.txt')
irg <- irg[!duplicated(irg$Symbol),]
idmapping <- read.table(file = "lianxi/PMID32798715/gencode.v22.annotation.gene.probeMap", 
                        sep = "\t", header = T, stringsAsFactors = F)
geneid <- data.frame(id = rownames(stad1), stringsAsFactors = F)
geneid2symbol <- dplyr::left_join(geneid, idmapping, by = "id")
stad1$symbol <- geneid2symbol$gene
stad1 <- stad1[which(stad1$symbol %in% irg$Symbol),]
stad1 <- stad1[order(stad1$symbol),]
stad2 <- aggregate(stad1[,c(1:(ncol(stad1)-1))], 
                   by = list(stad1$symbol), mean)
rownames(stad2) <- stad2[,1]
stad2 <- stad2[,-1]#
save(stad2,file = "lianxi/PMID32798715/STAD/stad_exp.RData")
```

## ESCA

```R
esca <- data.table::fread('lianxi/PMID32798715/ESCA/TCGA-ESCA.htseq_fpkm.tsv')
surv <- data.table::fread('lianxi/PMID32798715/ESCA/TCGA-ESCA.survival.tsv')
surv <- surv[,-3]
pheno <- data.table::fread('lianxi/PMID32798715/ESCA/TCGA-ESCA.GDC_phenotype.tsv')
pheno <- pheno[,c("submitter_id.samples","age_at_initial_pathologic_diagnosis","gender.demographic","pathologic_M",
                  "pathologic_N","pathologic_T","tumor_stage.diagnoses")]
#pheno <- na.omit(pheno)
pheno <- pheno[which(pheno$pathologic_T != "TX"),]
pheno <- pheno[which(pheno$pathologic_N != "NX"),]
table(pheno$pathologic_M)
pheno <- pheno[which(pheno$pathologic_M != ""),]
table(pheno$pathologic_N)
pheno <- pheno[which(pheno$pathologic_N != ""),]
table(pheno$pathologic_T)
pheno <- pheno[which(pheno$pathologic_T != "T0"),]
table(pheno$tumor_stage.diagnoses)
pheno <- pheno[which(pheno$tumor_stage.diagnoses != "not reported"),]
pheno[which(pheno$pathologic_N == "N3a"),5] <- "N3"
pheno[which(pheno$pathologic_N == "N3b"),5] <- "N3"
pheno[which(pheno$pathologic_T == "T1a"),6] <- "T1"
pheno[which(pheno$pathologic_T == "T1b"),6] <- "T1"
pheno[which(pheno$pathologic_T == "T2a"),6] <- "T2"
pheno[which(pheno$pathologic_T == "T2b"),6] <- "T2"
pheno[which(pheno$pathologic_M == "M1a"),4] <- "M1"
pheno[which(pheno$pathologic_M == "M1b"),4] <- "M1"
pheno[which(pheno$pathologic_N == "N1b"),5] <- "N1"
pheno[which(pheno$pathologic_N == "N1a"),5] <- "N1"
pheno[which(pheno$pathologic_N == "N1c"),5] <- "N1"
pheno[which(pheno$pathologic_N == "N2a"),5] <- "N2"
pheno[which(pheno$pathologic_N == "N2b"),5] <- "N2"
pheno[which(pheno$pathologic_T == "T4a"),6] <- "T4"
pheno[which(pheno$pathologic_T == "T4b"),6] <- "T4"

pheno[which(pheno$tumor_stage.diagnoses == "stage i"),7] <- "I"
pheno[which(pheno$tumor_stage.diagnoses == "stage ia"),7] <- "I"
pheno[which(pheno$tumor_stage.diagnoses == "stage ib"),7] <- "I"

pheno[which(pheno$tumor_stage.diagnoses == "stage ii"),7] <- "II"
pheno[which(pheno$tumor_stage.diagnoses == "stage iia"),7] <- "II"
pheno[which(pheno$tumor_stage.diagnoses == "stage iib"),7] <- "II"
pheno[which(pheno$tumor_stage.diagnoses == "stage iic"),7] <- "II"
pheno[which(pheno$tumor_stage.diagnoses == "stage iii"),7] <- "III"
pheno[which(pheno$tumor_stage.diagnoses == "stage iiia"),7] <- "III"
pheno[which(pheno$tumor_stage.diagnoses == "stage iiib"),7] <- "III"
pheno[which(pheno$tumor_stage.diagnoses == "stage iiic"),7] <- "III"
pheno[which(pheno$tumor_stage.diagnoses == "stage iv"),7] <- "IV"
pheno[which(pheno$tumor_stage.diagnoses == "stage iva"),7] <- "IV"
exp_sample <- colnames(esca)[-1]
surv_sample <- surv$sample
pheno_sample <- pheno$submitter_id.samples
jiaoji <- intersect(intersect(exp_sample,surv_sample),pheno_sample)

pheno1 <- pheno[which(pheno$submitter_id.samples %in% jiaoji),]
colnames(pheno1) <- c("sample_id","age","gender","M","N","T","stage")
esca <- tibble::column_to_rownames(esca,"Ensembl_ID")
esca1 <- esca[,jiaoji] 

surv1 <- surv[which(surv$sample %in% jiaoji),]

irg <- data.table::fread('D:/linshiwenjian/r-daima/免疫基因.txt')
irg <- irg[!duplicated(irg$Symbol),]
idmapping <- read.table(file = "lianxi/PMID32798715/gencode.v22.annotation.gene.probeMap", 
                        sep = "\t", header = T, stringsAsFactors = F)
geneid <- data.frame(id = rownames(esca1), stringsAsFactors = F)
geneid2symbol <- dplyr::left_join(geneid, idmapping, by = "id")
esca1$symbol <- geneid2symbol$gene
esca1 <- esca1[which(esca1$symbol %in% irg$Symbol),]
esca1 <- esca1[order(esca1$symbol),]
esca2 <- aggregate(esca1[,c(1:(ncol(esca1)-1))], 
                   by = list(esca1$symbol), mean)
rownames(esca2) <- esca2[,1]
esca2 <- esca2[,-1]#
save(esca2,file = "lianxi/PMID32798715/ESCA/esca_exp.RData")
```

## COAD

```R
coad <- data.table::fread('lianxi/PMID32798715/COAD/TCGA-COAD.htseq_fpkm.tsv')
surv <- data.table::fread('lianxi/PMID32798715/COAD/TCGA-COAD.survival.tsv')
surv <- surv[,-3]
pheno <- data.table::fread('lianxi/PMID32798715/COAD/TCGA-COAD.GDC_phenotype.tsv')
pheno <- pheno[,c("submitter_id.samples","age_at_initial_pathologic_diagnosis","gender.demographic","pathologic_M",
                  "pathologic_N","pathologic_T","tumor_stage.diagnoses")]
#pheno <- na.omit(pheno)
pheno <- pheno[which(pheno$pathologic_T != "TX"),]
pheno <- pheno[which(pheno$pathologic_N != "NX"),]
table(pheno$pathologic_M)
pheno <- pheno[which(pheno$pathologic_M != ""),]
table(pheno$pathologic_N)
pheno <- pheno[which(pheno$pathologic_N != ""),]
table(pheno$pathologic_T)
pheno <- pheno[which(pheno$pathologic_T != "Tis"),]
table(pheno$tumor_stage.diagnoses)
pheno <- pheno[which(pheno$tumor_stage.diagnoses != "not reported"),]
pheno[which(pheno$pathologic_N == "N3a"),5] <- "N3"
pheno[which(pheno$pathologic_N == "N3b"),5] <- "N3"
pheno[which(pheno$pathologic_T == "T1a"),6] <- "T1"
pheno[which(pheno$pathologic_T == "T1b"),6] <- "T1"
pheno[which(pheno$pathologic_T == "T2a"),6] <- "T2"
pheno[which(pheno$pathologic_T == "T2b"),6] <- "T2"
pheno[which(pheno$pathologic_M == "M1a"),4] <- "M1"
pheno[which(pheno$pathologic_M == "M1b"),4] <- "M1"
pheno[which(pheno$pathologic_N == "N1b"),5] <- "N1"
pheno[which(pheno$pathologic_N == "N1a"),5] <- "N1"
pheno[which(pheno$pathologic_N == "N1c"),5] <- "N1"
pheno[which(pheno$pathologic_N == "N2a"),5] <- "N2"
pheno[which(pheno$pathologic_N == "N2b"),5] <- "N2"
pheno[which(pheno$pathologic_T == "T4a"),6] <- "T4"
pheno[which(pheno$pathologic_T == "T4b"),6] <- "T4"

pheno[which(pheno$tumor_stage.diagnoses == "stage i"),7] <- "I"
pheno[which(pheno$tumor_stage.diagnoses == "stage ia"),7] <- "I"
pheno[which(pheno$tumor_stage.diagnoses == "stage ib"),7] <- "I"

pheno[which(pheno$tumor_stage.diagnoses == "stage ii"),7] <- "II"
pheno[which(pheno$tumor_stage.diagnoses == "stage iia"),7] <- "II"
pheno[which(pheno$tumor_stage.diagnoses == "stage iib"),7] <- "II"
pheno[which(pheno$tumor_stage.diagnoses == "stage iic"),7] <- "II"
pheno[which(pheno$tumor_stage.diagnoses == "stage iii"),7] <- "III"
pheno[which(pheno$tumor_stage.diagnoses == "stage iiia"),7] <- "III"
pheno[which(pheno$tumor_stage.diagnoses == "stage iiib"),7] <- "III"
pheno[which(pheno$tumor_stage.diagnoses == "stage iiic"),7] <- "III"
pheno[which(pheno$tumor_stage.diagnoses == "stage iv"),7] <- "IV"
pheno[which(pheno$tumor_stage.diagnoses == "stage iva"),7] <- "IV"
pheno[which(pheno$tumor_stage.diagnoses == "stage ivb"),7] <- "IV"
exp_sample <- colnames(coad)[-1]
surv_sample <- surv$sample
pheno_sample <- pheno$submitter_id.samples
jiaoji <- intersect(intersect(exp_sample,surv_sample),pheno_sample)

pheno1 <- pheno[which(pheno$submitter_id.samples %in% jiaoji),]
colnames(pheno1) <- c("sample_id","age","gender","M","N","T","stage")
coad <- tibble::column_to_rownames(coad,"Ensembl_ID")
coad1 <- coad[,jiaoji] 

surv1 <- surv[which(surv$sample %in% jiaoji),]

irg <- data.table::fread('D:/linshiwenjian/r-daima/免疫基因.txt')
irg <- irg[!duplicated(irg$Symbol),]
idmapping <- read.table(file = "lianxi/PMID32798715/gencode.v22.annotation.gene.probeMap", 
                        sep = "\t", header = T, stringsAsFactors = F)
geneid <- data.frame(id = rownames(coad1), stringsAsFactors = F)
geneid2symbol <- dplyr::left_join(geneid, idmapping, by = "id")
coad1$symbol <- geneid2symbol$gene
coad1 <- coad1[which(coad1$symbol %in% irg$Symbol),]
coad1 <- coad1[order(coad1$symbol),]
coad2 <- aggregate(coad1[,c(1:(ncol(coad1)-1))], 
                   by = list(coad1$symbol), mean)
rownames(coad2) <- coad2[,1]
coad2 <- coad2[,-1]#
save(coad2,file = "lianxi/PMID32798715/COAD/coad_exp.RData")
```

