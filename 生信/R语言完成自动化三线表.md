#'http://bio-info-trainee.com/tmp/TCGA-LUAD-phe_clinical_tables.Rdata'
load('lianxi/TCGA-LUAD-phe_clinical_tables.Rdata')
clinical_info=phe
head(clinical_info) 
#首先对需要观测的临床特质值进行重新编码 
clinical_info$age<-as.numeric(clinical_info$age)
clinical_info$AGE<-factor(ifelse(clinical_info$age>60,'>60','<=60'),ordered = T)
clinical_info$gender<-factor(toupper(clinical_info$gender),levels=c("MALE", "FEMALE"),ordered = T)
clinical_info$stage<-factor(toupper(clinical_info$stage),ordered = T)
clinical_info$t<-factor(clinical_info$t,ordered = T)
clinical_info$n<-factor(clinical_info$n,ordered = T)
clinical_info$m<-factor(clinical_info$m,ordered = T) 
clinical_info$vital_status<-factor(toupper(clinical_info$vital_status),ordered = T)
clinical_info$race<-factor(clinical_info$race,ordered = T)

#去除不需要的临床信息

clinical_info=clinical_info[,-c(1,10:12)]
dput(names(clinical_info))

#Vector of variables to summarize

myVars <- dput(names(clinical_info))

#Vector of categorical variables that need transformation

catVars <- myVars[c(1,2,4:8,10)]

##三线表类型之一  切割数据 
library(caret)
set.seed(123456789)
sam<- createDataPartition(clinical_info$vital_status, p = .5,list = FALSE)
train <- clinical_info[sam,]
test <- clinical_info[-sam,]
#查看两组一些临床参数切割比例
prop.table(table(train$stage))
prop.table(table(test$stage))
#添加分组
train$group<-'training datasets'
test$group<-'testing datasets'
clinical_info<-rbind(train,test)
clinical_info$group<-factor(clinical_info$group)
##生成三线表
vars <-colnames(clinical_info)[c(2:9,12,14,15)]
library(tableone)

#最重要的三线表通常是以训练集和数据集来区分：group

tb_group<-CreateTableOne(vars = myVars, strata = c("group"), data = clinical_info,
                         factorVars = catVars) 
tab1<-print(tb_group, nonnormal = c('age','time'),
            exact = c(myVars,'AGE'), smd = TRUE)
summary(tab1)
tab_out<-print(tb_group, catDigits = 1, contDigits = 2, pDigits = 3,
               quote = FALSE, missing = T, explain = TRUE, printToggle = TRUE,
               test = TRUE, smd = T, noSpaces = FALSE, padColnames = FALSE,
               varLabels = FALSE, format = c("fp", "f", "p", "pf")[1],
               showAllLevels = FALSE, cramVars = NULL, dropEqual = FALSE,
               exact = NULL, nonnormal = NULL, minMax = FALSE)

write.csv(tab_out, file = "lianxi/TCGA-LUAD-phe_clinical_tables1.csv")