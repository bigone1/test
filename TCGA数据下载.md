***

title：TCGA数据下载

author：周林

date：2020.07.02

***

### 一、GDC下载：https://portal.gdc.cancer.gov/repository

1. 下载表达数据和临床数据的manifest文件以及metadata文件：

   ![Error](https://github.com/bigone1/test/blob/master/Screenshots/9.png)
   
2. 下载数据：

   ```shell
   /home/zhoulin/gdc-client download -m gdc_manifest.LIHC_mirna.txt -d mirna/
   
   /home/zhoulin/gdc-client download -m gdc_manifest.LIHC_clinical.txt -d clinical/
   ```

3. 表达矩阵和临床信息整理：

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
   mis <- dir('mirna/',pattern = '*tification.txt$',recursive = T)
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

   