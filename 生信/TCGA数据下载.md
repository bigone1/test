## 一、GDC下载：https://portal.gdc.cancer.gov/repository

### 下载表达数据和临床数据的manifest文件以及metadata文件

![Error](https://github.com/bigone1/test/blob/master/Screenshots/9.png)

### 下载数据

```shell
/home/zhoulin/gdc-client download -m gdc_manifest.LIHC_mirna.txt -d mirna/

/home/zhoulin/gdc-client download -m gdc_manifest.LIHC_clinical.txt -d clinical/
```

### 表达矩阵和临床信息整理

数据位于：https://github.com/bigone1/test/blob/master/data/LIHC

```R
library(XML)

#整理临床信息
#单个样本的临床信息#####
result <- xmlParse('clinical/004d6594-95ce-494c-9760-828f9885bbae/nationwidechildrens.org_clinical.TCGA-DD-A1EA.xml')
rootnode <- xmlRoot(result)
rootsize <- xmlSize(rootnode)#两个节点
print(rootnode[2])#病人信息在第二部分
xmldataframe <- xmlToDataFrame(rootnode[2])
#####
#处理所有样本的临床信息#####
xmls <- dir('clinical/',pattern = '*.xml$',recursive = T)
td <- function(x){
  result <- xmlParse(file.path('clinical/',x))
  rootnode <- xmlRoot(result)
  xmldataframe <- xmlToDataFrame(rootnode[2])[c(
    'bcr_patient_barcode',
    'vital_status',
    'days_to_death',
    'days_to_last_followup',
    'race_list',
    'days_to_birth',
    'gender',
    'stage_event'
  )]
  return(t(xmldataframe))
}
cl <- lapply(xmls, td)
cl_df <- t(do.call(cbind,cl))
clinical <- data.frame(cl_df)
rownames(clinical) <- stringr::str_to_upper(clinical$bcr_patient_barcode)
clinical <- clinical[,-1]
#####


#整理表达矩阵
mis <- dir('mirna/',pattern = '*tification.txt$',recursive = T)#注意文件名
ex <- function(x){
  result <- read.table(file.path('mirna/',x),sep='\t',header = T)[,1:2]
  return(result)
}
mi <- lapply(mis, ex)
mi_df <- t(do.call(cbind,mi))
colnames(mi_df) <-mi_df[1,]
mi_df <- mi_df[seq(2,nrow(mi_df),2),]
mi_df <- apply(mi_df,2,as.numeric)
#用metadata文件给表达矩阵加上样本名
meta <- jsonlite::fromJSON('metadata.LIHC.json')
entity <- meta$associated_entities
jh <- function(x){
  as.character(x[1])
}
ID <- sapply(entity,jh)
options(stringsAsFactors = F)
file2id <- data.frame(file_name = meta$file_name,ID = ID)
mis2 <- stringr::str_split(mis,"/",simplify = T)[,2]
row_tcga <- file2id[match(mis2,file2id$file_name),]
rownames(mi_df) <- row_tcga$ID
expr <- t(mi_df) 
expr <- expr[apply(expr,1,function(x){
  sum(x >1) > 9
}),]
group_list <- ifelse(as.numeric(substr(colnames(expr),14,15)) < 10, 'tumor','normal')
group_list <- factor(group_list,levels = c('normal','tumor'))
save(expr,clinical,group_list,file="gdc.Rdata")
```

# 二、使用TCGAbiolinks批量下载TCGA的表达量数据

```R
#使用TCGAbiolinks批量下载TCGA的表达量数据
library(TCGAbiolinks)
#查看可以下载的项目
project <- getGDCprojects()
#提取TCGA项目
library(dplyr)
projects <- project %>% 
  as.data.frame() %>% 
  select(project_id,tumor) %>% 
  filter(grepl(pattern="TCGA",project_id))
#下载表达矩阵，分四步
#一：使用GDCquery函数来查询下载信息
query.exp <- GDCquery(project = "TCGA-OV", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts")
#二：使用GDCdownload函数来下载数据
GDCdownload(query.exp)
#三：使用GDCprepare批量读取数据，同时保存为Rdata数据格式
pre.exp <- GDCprepare(query = query.exp, save = TRUE, 
                      save.filename = "TCGA_OV_RNAseq_counts.Rdata")
#四：使用assay函数提取表达量信息
countsdata <- SummarizedExperiment::assay(pre.exp)


#批量下载TCGA的数据
dir.create("TCGA_RNA_data")
for (i in 1:nrow(projects)) {
  ## 0.运行信息
  print(paste0("Downloading number ",i,",project name: ",projects$project_id[i]))
  ## 1.查询信息
  query.exp = GDCquery(project = projects$project_id[i], 
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "HTSeq - Counts")
  ## 2.正式下载
  GDCdownload(query.exp)
  ## 3.多个数据合并
  pre.exp = GDCprepare(query = query.exp)
  ## 4.提取表达量数据
  countsdata = SummarizedExperiment::assay(pre.exp)
  ## 5.保存数据
  save(countsdata,file = paste0("TCGA_RNA_data/",projects$project_id[i],".Rdata"))
}
```

# 官网下载

### 一：TCGA数据下载：[https://portal.gdc.cancer.gov](https://portal.gdc.cancer.gov/)，下载mnifest,samplesheet文件，通过以下命令下载数据

```powershell
gdc-client download -m gdc_manifest.txt
```

### 二：gdc_sample_sheet.tsv文件删选样本

```R
samsheet <- read.csv('gdc_sample_sheet.2019-09-29.tsv',sep='\t',header = T)
samsheet_filter<-samsheet[substr(samsheet$Sample.ID,14,17)=='01A',]
library(dplyr)
c<-filter(samsheet_filter,!duplicated(samsheet_filter$Sample.ID))
write.csv(c,'sample-sheet.csv',row.names = F)
```

### 三：使用samplesheet文件根据Sample Type删除正常组织样本，提取File ID、File Name、Sample ID三个信息存为txt文件：sample-sheet.txt

### 四：创建file文件夹，用下面的脚本extract.sh，执行：bash extract.sh sample-sheet.txt，获得count文件，删除不正确的count文件

```perl
 #!/bin/bash
 ?
 cat $1 |while read line
 do
  arr=($line)
  filename=${arr[1]}
  folder=${arr[0]}
  submitterid=${arr[2]}
 gunzip -c ./${folder}/${filename} > ./file/${submitterid:0:-1}.count
 done
```

### 五：下面的脚本合并count文件，执行：perl merge.pl ./file/ > merge.count

```perl
#!/bin/perl -w
use strict;

my $path = shift @ARGV;

opendir DIR, $path or die;
my @dir = readdir DIR;

my $header;
my @sample;
my %hash;
my %rm_rep;

foreach my $file (@dir) {
    if ($file =~ /^(\w+.*-\d+)[A-Z]\.count/) {
        my $sample_id = $1;
        next unless ($sample_id =~ /01$/ || $sample_id =~ /11$/);
        if ($rm_rep{$sample_id}) {
            next;
        }else{
            $rm_rep{$sample_id} = 1;
        }
        push @sample, $sample_id;
        $header .= "\t$sample_id";

        open my $fh, $path.$file or die;
        while (<$fh>) {
            chomp;
            next if ($_ =~ /^_/);
            my @array = split /\t/, $_;
            $hash{$array[0]} -> {$sample_id} = $array[1];
        }
        close $fh;
    }
}

print "$header\n";
map{
    my $gene = $_;
    print "$gene";
    foreach my $file (@sample) {
        print "\t".$hash{$gene} -> {$file};
    }
    print "\n";
}keys %hash;
```

### 六：以下脚本：zhushi.R 对合并的矩阵注释，获得GeneSymbol文件merge.txt

```R
raw_count <- read.table('merge.count',sep = '\t',header = T)
#raw_count_filt <- raw_count[-(nrow(raw_count)-4):-nrow(raw_count),]
ENSEMBL <- gsub("\\.\\d*", "", raw_count$X)
row.names(raw_count) <- ENSEMBL
raw_count_filt <- cbind(ENSEMBL,raw_count)
colnames(raw_count_filt)[1] <- c("ensembl_gene_id")
library('biomaRt')
library("curl")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id <- row.names(raw_count_filt)
options(timeout = 4000000)
hg_symbols<- getBM(attributes=c('ensembl_gene_id','hgnc_symbol',"chromosome_name", "start_position","end_position", "band"), filters= 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
readcount <- merge(raw_count_filt, hg_symbols, by="ensembl_gene_id")
readcount <- readcount[readcount$hgnc_symbol!="",]
readcount <- readcount[,-2]
readcount<-readcount[,-(ncol(readcount)-3):-ncol(readcount)]
readcount$ensembl_gene_id=readcount$hgnc_symbol
readcount<-readcount[,-ncol(readcount)]
write.table(readcount,'merge.txt',sep = '\t',row.names = F)
```

# 官网下载整合-差异分析

## 一、下载cart和Metadata数据

![Error](https://github.com/bigone1/test/blob/master/Screenshots/19.png)

## 二、解压cart文件到单独文件夹

```R
###################################数 据 处 理#####################################
options(stringsAsFactors = F)
dir.create('samplefiles')
filepath <- dir('./mrna',full.names = TRUE)
for(wd in filepath){
  files <- dir(path = wd,pattern = "gz$")#查看满足条件文件
  fromfilepath <- paste(wd,'/',files,sep="")
  tofilepath <- paste("./samplefiles/",files,sep="")
  file.copy(fromfilepath,tofilepath)
}
################################### 解压所有文件并删除原文件 ########################
setwd('./samplefiles/')
countsFiles <- dir(path = './',pattern = 'gz$')#查看满足条件文件
library(R.utils)
sapply(countsFiles,gunzip)#解压函数gunzip需要R.utils包

################################### 处理json文件 ########################
library(rjson)
metadata_json_File <- fromJSON(file = 'D:/linshiwenjian/r-daima/linshi/luad/metadata.LUAD.json')
json_File_Info <- data.frame(fileName=c(),TCGA_Barcode=c())
for(i in 1:length(metadata_json_File)){
  TCGA_Barcode <- metadata_json_File[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
  file_name <- metadata_json_File[[i]][["file_name"]]
  json_File_Info <- rbind(json_File_Info,data.frame(filesName=file_name,TCGA_Barcode=TCGA_Barcode))
}
rownames(json_File_Info) <- json_File_Info[,1]
write.csv(json_File_Info,file = "json_File_Info.csv")


################################### 获取Counts矩阵 ########################
filesName_To_TCGA_BarcodeFile <- json_File_Info[-1]
countsFileNames <- dir(pattern = "txt$")
allSampleRawCounts <- data.frame()
for(txtFile in countsFileNames){
  #每个循环读取一个文件
  SampleCounts <- read.table(txtFile,header = FALSE)
  rownames(SampleCounts) <- SampleCounts[,1]
  SampleCounts <- SampleCounts[-1]
  #根据filesName_To_TCGA_BarcodeFile文件中的文件名称与barcode对应关系，命名列名
  colnames(SampleCounts) <- filesName_To_TCGA_BarcodeFile[paste(txtFile,".gz",sep=""),1]#
  if(dim(allSampleRawCounts)[1]==0){
    allSampleRawCounts <- SampleCounts
  }
  else{
    allSampleRawCounts <- cbind(allSampleRawCounts,SampleCounts)
  }
}
write.csv(allSampleRawCounts,file='../allSampleRawCounts.csv')
ensemble_id <- substr(row.names(allSampleRawCounts),1,15)
rownames(allSampleRawCounts) <- ensemble_id
#RawCounts.csv与allSampleRawCounts.csv文件的区别是行名的ensembl去掉了版本号
write.csv(allSampleRawCounts,file='../RawCounts.csv')

################################### ID转换 ########################
#在GENCODE上下载注释文件gencode.v33lift37.annotation.gtf
RawCounts <- allSampleRawCounts
Ensembl_ID <- data.frame(Ensembl_ID=row.names(RawCounts))
rownames(Ensembl_ID) <- Ensembl_ID[,1]
RawCounts <- cbind(Ensembl_ID,RawCounts)
#该函数，通过gtf文件获取Ensembl_ID与基因名称的对应关系
get_map <- function(input){
  if(is.character(input)){
    if(!file.exists(input)){
      stop("Bad input file.")
    }
    message("Treat input as file")
    input <- data.table::fread(input,header = FALSE)
  }else{
    data.table::setDT(input)
  }
  input <- input[input[[3]]=='gene',]
  pattern_id <- ".*gene_id \"([^;]+)\";.*"
  pattern_name <- ".*gene_name \"([^;]+)\";.*"
  gene_id <- sub(pattern_id,"\\1",input[[9]])
  gene_name <- sub(pattern_name,"\\1",input[[9]])
  Ensembl_ID_To_Genename <- data.frame(gene_id=gene_id,
                                       gene_name=gene_name,
                                       stringsAsFactors = FALSE)
  return(Ensembl_ID_To_Genename)
}

Ensembl_ID_To_Genename <- get_map('../gencode.v33lift37.annotation.gtf')
gtf_Ensembl_ID <- substr(Ensembl_ID_To_Genename[,1],1,15)
Ensembl_ID_To_Genename <- data.frame(gtf_Ensembl_ID,Ensembl_ID_To_Genename[,2])
colnames(Ensembl_ID_To_Genename) <- c("Ensembl_ID","gene_id")
write.csv(Ensembl_ID_To_Genename,file='../Ensembl_ID_To_Genename.csv')
#融合数据
mergeRawCounts <- merge(Ensembl_ID_To_Genename,RawCounts,by='Ensembl_ID')
#按照gene_id列进行排序
mergeRawCounts <- mergeRawCounts[order(mergeRawCounts[,"gene_id"]),]
#根据gene_id列建立索引
index <- duplicated(mergeRawCounts$gene_id)
#提取不重复的基因
mergeRawCounts <- mergeRawCounts[!index,]
rownames(mergeRawCounts) <- mergeRawCounts[,"gene_id"]
LUAD_Cunts_expMatrix <- mergeRawCounts[,-c(1:2)]
write.csv(LUAD_Cunts_expMatrix,file='../LUAD_Cunts_expMatrix.csv')

################################### 差异分析 ########################
#读取数据：Counts <- read.csv('../LUAD_Cunts_expMatrix.csv',header=T,row.names=1)
library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
#从query中获取结果表，可以是带有cols参数的列，也有rows参数返回若干行
samplesDown <- getResults(query,cols = c("cases"))
#肿瘤样本的barcode
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,typesample = "TP")
#正常样本的barcode
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,typesample = "NT")
#重新排序，正常样本在前
Counts <- data.frame(c(LUAD_Cunts_expMatrix[,dataSmNT],LUAD_Cunts_expMatrix[,dataSmTP]))
rownames(Counts) <- row.names(LUAD_Cunts_expMatrix)
colnames(Counts) <- c(dataSmNT,dataSmTP)

#edgeR中，1代表control样本，2代表case样本
#######################方法一：edgeR
library(edgeR)
group <- c(rep(1,59),rep(2,533))
y <- DGEList(counts = Counts,group = group)
#数据过滤
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
#计算标准化因子
y <- calcNormFactors(y)
#计算离散度
y <- estimateDisp(y)
#显著性检验
et <- exactTest(y)

et <- topTags(et,n=100000)
et <- as.data.frame(et)
et <- cbind(rownames(et),et)
colnames(et) <- c("gene_id","log2FoldChange","log2CPM","Pvalue","FDR")
write.table(et,file = '../all_LUAD_DEG.xls',sep='\t',col.names = TRUE,
            row.names = FALSE,quote = FALSE,na="")
etSig <- et[which(et$Pvalue < 0.05 & abs(et$log2FoldChange) > 1),]
etSig[which(etSig$log2FoldChange > 0),"up_down"] <- "up"
etSig[which(etSig$log2FoldChange < 0),"up_down"] <- "down"
write.table(etSig,file = '../LUAD_DEG.xls',sep='\t',col.names = TRUE,
            row.names = FALSE,quote = FALSE,na="")

#############方法二：DESeq2
library(DESeq2)
DESeq2group <- c(rep("control",59),rep("case",533))
colData <- data.frame(row.names = names(Counts),
                      condition=factor(DESeq2group,levels = c("control","case")))
dds <- DESeqDataSetFromMatrix(countData = Counts,colData = colData,
                              design = ~ condition)
dds <- DESeq(dds)
sizeFactors(dds)
res <- results(dds)
res <- as.data.frame(res)
res <- cbind(rownames(res),res)
colnames(res) <- c("gene_id","baseMean","log2FoldChange","lfcSE","stat",
                   "pval","padj")
write.table(res,file = '../case-vs-control-all-DESeq2.gene.xls',sep='\t',col.names = TRUE,
            row.names = FALSE,quote = FALSE,na="")
resSig <- res[which(res$pval < 0.05 & abs(res$log2FoldChange) > 1),]
resSig[which(resSig$log2FoldChange > 0),"up_down"] <- "up"
resSig[which(resSig$log2FoldChange < 0),"up_down"] <- "down"
write.table(resSig,file = '../case-vs-control-diff-pval-0.05-FC-2-DESeq2.gene.xls',sep='\t',col.names = TRUE,
            row.names = FALSE,quote = FALSE,na="")
```

# 临床数据下载整理

```R
rm(list = ls())
options(stringsAsFactors = F)
project <- 'TCGA-COAD'
clinicaldir <- "D:/download/gdc_download_20201111_062248.922899/"

library(XML)
library(methods)
xmlFileName <- dir(path = clinicaldir,full.names = T,
                   pattern = "xml$",recursive = T)
AllPaticliniList <- lapply(xmlFileName,
                           function(x){
                             result <- xmlParse(file = x)
                             rootNode <- xmlRoot(result)
                             xmldataframe <- xmlToDataFrame(rootNode[2])
                             return(t(xmldataframe))
                           })
AllPaticliniList <- t(do.call(cbind,AllPaticliniList))
pclinData <- data.frame(Barcode=AllPaticliniList[,"bcr_patient_barcode"],
                        Vital_status=AllPaticliniList[,"vital_status"],
                        days_to_death=AllPaticliniList[,"days_to_death"],
                        days_to_last_known_alive=AllPaticliniList[,"days_to_last_known_alive"],
                        lastFollowupTime=AllPaticliniList[,"days_to_last_followup"],
                        Stage=AllPaticliniList[,"stage_event"])
rownames(pclinData) <- pclinData[,"Barcode"]
splitStage <- function(Stage){
  stage <- unlist(strsplit(Stage,split = "[S]"))[2]
  stage <- unlist(strsplit(stage,split = "T"))[1]
  stage <- paste("S",stage,sep = "")
  
  Tstage <- unlist(strsplit(Stage,split = "[T]"))[2]
  Tstage <- unlist(strsplit(Tstage,split = "N"))[1]
  Tstage <- paste("T",Tstage,sep = "")
  
  Nstage <- unlist(strsplit(Stage,split = "[N]"))[2]
  Nstage <- unlist(strsplit(Nstage,split = "M"))[1]
  Nstage <- paste("N",Nstage,sep="")
  
  Mstage <- unlist(strsplit(Stage,split = "[M]"))[2]
  Mstage <- paste("M",Mstage,sep="")
  
  return(data.frame(stage,Tstage,Nstage,Mstage))
}
SurvivalTime <- c()
Stage <- c()
T_stage <- c()
N_stage <- c()
M_stage <- c()
for(id in rownames(pclinData)){
  days_to_death = as.numeric(pclinData[id,"days_to_death"])
  lastFollowupTime=as.numeric(pclinData[id,"lastFollowupTime"])
  Patientsurvtime <- c()
  if(!is.na(days_to_death)){
    Patientsurvtime <- c(days_to_death)
  }else{Patientsurvtime <- c(lastFollowupTime)}
  SurvivalTime <- c(SurvivalTime,Patientsurvtime)
  TNMStage <- splitStage(Stage = pclinData[id,"Stage"])
  Stage <- c(Stage,TNMStage[1,"stage"])
  T_stage <- c(T_stage,TNMStage[1,"Tstage"])
  N_stage <- c(N_stage,TNMStage[1,"Nstage"])
  M_stage <- c(M_stage,TNMStage[1,"Mstage"])
}
getClinicalData <- function(pclinData){
  tempdata <- pclinData
  finalpclinical <- data.frame(Barcode=tempdata[,"Barcode"],
                               vital_status=tempdata[,"Vital_status"],
                               SurvivalTime=SurvivalTime,
                               AllStage=tempdata[,"Stage"],
                               Stage=Stage,
                               T_stage=T_stage,
                               N_stage=N_stage,
                               M_stage=M_stage)
  return(finalpclinical)
}
tidyclidata <- getClinicalData(pclinData)

library(GDCRNATools)
gdcClinicalDownload(project.id = project,
                    write.manifest = FALSE,
                    method = 'gdc-client',
                    directory = 'D:/File/Rfile/clinical')
```

# ID转换

```R
library(biomaRt)
listMarts()#查看支持的数据库
ensemble=useMart("ensembl")#选择一个数据库

#选择数据集
ensemble=useDataset('hsapiens_gene_ensembl',mart = ensemble)
#我们需要使用的属性有NM,geneSymbol
NMID=c("NM_001193582","NM_004951","NM_003489")
query=getBM(attributes = c("hgnc_symbol",
"refseq_mrna"),
            filters = ('refseq_mrna'),
            values = list(NMID),
            mart = ensemble)

#query
#  hgnc_symbol  refseq_mrna
#1       NRCAM NM_001193582
#2       NRIP1    NM_003489
#3      GPR183    NM_004951

gse76427 <- getGEO("GSE76427",AnnotGPL=F,getGPL=F)
gse76427_exp <- exprs(gse76427[[1]])
library(biomaRt)
ensem <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl",mart=ensem)
query76427 <- getBM(attributes=c("illumina_humanht_12_v4","entrezgene_id"),filters=("illumina_humanht_12_v4"),values=rownames(gse76427_exp),mart=mart)
prob76427 <- data.frame(illumina_humanht_12_v4=rownames(gse76427_exp),stringsAsFactors=F)
prob76427s <- dplyr::left_join(prob76427,query76427,by='illumina_humanht_12_v4')
prob76427sn <- na.omit(prob76427s)
j76427 <- intersect(prob76427sn$illumina_humanht_12_v4,rownames(gse76427_exp))
gse76427_exp <- gse76427_exp[j76427,]
gse76427_exp <- as.data.frame(gse76427_exp)
prob76427snc <-dplyr::filter(prob76427sn,!duplicated(prob76427sn$illumina_humanht_12_v4))
gse76427_exp$entr <- prob76427snc$entrezgene_id
gse76427_exp <- gse76427_exp[order(gse76427_exp$entr),]
gse76427_exp1 <- aggregate(gse76427_exp[,c(1:(ncol(gse76427_exp)-1))],by=list(gse76427_exp$entr),mean)
write.table(gse76427_exp1,file='gse76427.txt',sep='\t')

gse36376 <- getGEO("GSE36376",AnnotGPL=F,getGPL=F)
gse36376_exp <- exprs(gse36376[[1]])
library(biomaRt)
ensem <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl",mart=ensem)
query36376 <- getBM(attributes=c("illumina_humanht_12_v4","entrezgene_id"),filters=("illumina_humanht_12_v4"),values=rownames(gse36376_exp),mart=mart)
prob36376 <- data.frame(illumina_humanht_12_v4=rownames(gse36376_exp),stringsAsFactors=F)
prob36376s <- dplyr::left_join(prob36376,query36376,by='illumina_humanht_12_v4')
prob36376sn <- na.omit(prob36376s)
j36376 <- intersect(prob36376sn$illumina_humanht_12_v4,rownames(gse36376_exp))
gse36376_exp <- gse36376_exp[j36376,]
gse36376_exp <- as.data.frame(gse36376_exp)
prob36376snc <-dplyr::filter(prob36376sn,!duplicated(prob36376sn$illumina_humanht_12_v4))
gse36376_exp$entr <- prob36376snc$entrezgene_id
gse36376_exp <- gse36376_exp[order(gse36376_exp$entr),]
gse36376_exp1 <- aggregate(gse36376_exp[,c(1:(ncol(gse36376_exp)-1))],by=list(gse36376_exp$entr),mean)
write.table(gse36376_exp1,file='gse36376.txt',sep='\t')

gse84005 <- getGEO("GSE84005",AnnotGPL=F,getGPL=F)
gse84005_exp <- exprs(gse84005[[1]])
library(biomaRt)
ensem <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl",mart=ensem)
query84005 <- getBM(attributes=c("affy_huex_1_0_st_v2","entrezgene_id"),filters=("affy_huex_1_0_st_v2"),values=rownames(gse84005_exp),mart=mart)
prob84005 <- data.frame(illumina_humanht_12_v4=rownames(gse84005_exp),stringsAsFactors=F)
prob84005s <- dplyr::left_join(prob84005,query84005,by='illumina_humanht_12_v4')
prob84005sn <- na.omit(prob84005s)
j84005 <- intersect(prob84005sn$illumina_humanht_12_v4,rownames(gse84005_exp))
gse84005_exp <- gse84005_exp[j84005,]
gse84005_exp <- as.data.frame(gse84005_exp)
prob84005snc <-dplyr::filter(prob84005sn,!duplicated(prob84005sn$illumina_humanht_12_v4))
gse84005_exp$entr <- prob84005snc$entrezgene_id
gse84005_exp <- gse84005_exp[order(gse84005_exp$entr),]
gse84005_exp1 <- aggregate(gse84005_exp[,c(1:(ncol(gse84005_exp)-1))],by=list(gse84005_exp$entr),mean)
write.table(gse84005_exp1,file='gse84005.txt',sep='\t')

gse101685 <- getGEO("GSE101685",AnnotGPL=F,getGPL=F)
gse101685_exp <- exprs(gse101685[[1]])
library(biomaRt)
ensem <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl",mart=ensem)
query101685 <- getBM(attributes=c("affy_hg_u133_plus_2","entrezgene_id"),filters=("affy_hg_u133_plus_2"),values=rownames(gse101685_exp),mart=mart)
prob101685 <- data.frame(affy_hg_u133_plus_2=rownames(gse101685_exp),stringsAsFactors=F)
prob101685s <- dplyr::left_join(prob101685,query101685,by='affy_hg_u133_plus_2')
prob101685sn <- na.omit(prob101685s)
j101685 <- intersect(prob101685sn$affy_hg_u133_plus_2,rownames(gse101685_exp))
gse101685_exp <- gse101685_exp[j101685,]
gse101685_exp <- as.data.frame(gse101685_exp)
prob101685snc <-dplyr::filter(prob101685sn,!duplicated(prob101685sn$affy_hg_u133_plus_2))
gse101685_exp$entr <- prob101685snc$entrezgene_id
gse101685_exp <- gse101685_exp[order(gse101685_exp$entr),]
gse101685_exp1 <- aggregate(gse101685_exp[,c(1:(ncol(gse101685_exp)-1))],by=list(gse101685_exp$entr),mean)
write.table(gse101685_exp1,file='gse101685.txt',sep='\t')
```

illumina_humanht_12_v4

affy_huex_1_0_st_v2

affy_hg_u133_plus_2

entrezgene_id

ensembl_gene_id



```R
s012682 <- recount_url[which(recount_url$project == 'SRP012682'),]
metadata <- all_metadata('gtex')
ensg <- unlist(strsplit(rownames(rpkm),"[.]"))
ens <- ensg[seq(1,length(ensg),2)]
en <- data.frame(id=ens)
rpkm <- cbind(en,rpkm)
rpkm <- dplyr::filter(rpkm,!duplicated(rpkm$id))
query <- getBM(attributes=c("ensembl_gene_id","entrezgene_id"),filters=("ensembl_gene_id"),values=rownames(rpkm),mart=mart)
probe <- data.frame(ensembl_gene_id=rownames(rpkm),stringsAsFactors=F)
probe2id <- dplyr::left_join(probe,query,by='ensembl_gene_id')
probe2id <- dplyr::filter(probe2id,!duplicated(probe2id$ensembl_gene_id))
probe2id <- tibble::column_to_rownames(probe2id,"ensembl_gene_id")
rpkm <- rpkm[order(rpkm$entre),]
rpkm1 <- aggregate(rpkm[,c(1:(ncol(rpkm)-1))],by=list(rpkm$entre),mean)
write.table(rpkm1,file='gtex_liver.txt',sep='\t')
savehistory('his.txt')

library('org.Hs.eg.db')
gencode <- gsub('\\..*', '', names(recount_genes))
gene_info <- select(org.Hs.eg.db, gencode, c('ENTREZID', 'GENENAME', 'SYMBOL',
    'ENSEMBL'), 'ENSEMBL')

```

