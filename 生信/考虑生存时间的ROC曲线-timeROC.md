### 1.准备输入数据

输入数据是TCGA的表达矩阵expr、临床信息meta和group_list。保存为forest.Rdata了，在生信星球公众号后台聊天窗口回复“森林”即可获得。（为什么是森林呢，因为和随机森林图用的同一个数据，我又又又又懒得改了。）

```R
load("forest.Rdata")
exprSet = expr[,group_list=="tumor"]

## 必须保证生存资料和表达矩阵，两者一致
all(substring(colnames(exprSet),1,12)==meta$ID)
```

### 2.构建lasso回归模型

输入数据是表达矩阵(仅含tumor样本)和对应的生死。

```R
x=t(log2(exprSet+1))
y=meta$event
library(glmnet)
cv_fit <- cv.glmnet(x=x, y=y, nlambda = 1000,alpha = 1)
model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
length(choose_gene_min)
```

### 3.模型预测和评估

#### 3.1自己预测自己

```R
lasso.prob <- predict(cv_fit, newx=x , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
re=cbind(y ,lasso.prob)
new_dat=meta
library(timeROC)
library(survival)
library(survminer)
new_dat$fp=as.numeric(lasso.prob[,1])
with(new_dat,
     ROC <<- timeROC(T=time,#结局时间 
                     delta=event,#生存结局 
                     marker=fp,#预测变量 
                     cause=1,#阳性结局赋值，比如死亡与否
                     weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
                     times=c(60,100),#时间点，选取5年(60个月)和8年生存率
                     ROC = TRUE,
                     iid = TRUE)
)
auc_60 = ROC$AUC[[1]]
auc_100 = ROC$AUC[[2]]
dat = data.frame(tpr_60 = ROC$TP[,1],
                 fpr_60 = ROC$FP[,1],
                 tpr_100 = ROC$TP[,2],
                 fpr_100 = ROC$FP[,2])
library(ggplot2)

ggplot() + 
  geom_line(data = dat,aes(x = fpr_60, y = tpr_60),color = "blue") + 
  geom_line(data = dat,aes(x = fpr_100, y = tpr_100),color = "red")+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  annotate("text",x = .75, y = .25,
           label = paste("AUC of 60 = ",round(auc_60,2)),color = "blue")+
  annotate("text",x = .75, y = .15,label = paste("AUC of 100 = ",round(auc_100,2)),color = "red")+
  scale_x_continuous(name  = "fpr")+
  scale_y_continuous(name = "tpr")
```

### 4.切割数据构建模型并预测

#### 4.1 切割数据

用R包caret切割数据，生成的结果是一组代表列数的数字，用这些数字来给表达矩阵和meta取子集即可。

```R
library(caret)
set.seed(12345679)
sam<- createDataPartition(meta$event, p = .5,list = FALSE)
train <- exprSet[,sam]
test <- exprSet[,-sam]
train_meta <- meta[sam,]
test_meta <- meta[-sam,]
```

#### 4.2 切割后的train数据集建模

和上面的建模方法一样。

```R
#计算lambda
x = t(log2(train+1))
y = train_meta$event
cv_fit <- cv.glmnet(x=x, y=y, nlambda = 1000,alpha = 1)
#构建模型
model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
#挑出基因
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
length(choose_gene_min)
lasso.prob <- predict(cv_fit, newx=t(log2(test+1)), s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
re=cbind(test_meta$event ,lasso.prob)
new_dat = test_meta
library(timeROC)
library(survival)
library(survminer)
new_dat$fp=as.numeric(lasso.prob[,1])
with(new_dat,
     ROC <<- timeROC(T=time,#结局时间 
                     delta=event,#生存结局 
                     marker=fp,#预测变量 
                     cause=1,#阳性结局赋值，比如死亡与否
                     weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
                     times=c(60,100),#时间点，选取5年(60个月)和8年生存率
                     ROC = TRUE,
                     iid = TRUE)
)
auc_60 = ROC$AUC[[1]]
auc_100 = ROC$AUC[[2]]
dat = data.frame(tpr_60 = ROC$TP[,1],
                 fpr_60 = ROC$FP[,1],
                 tpr_100 = ROC$TP[,2],
                 fpr_100 = ROC$FP[,2])
library(ggplot2)

ggplot() + 
  geom_line(data = dat,aes(x = fpr_60, y = tpr_60),color = "blue") + 
  geom_line(data = dat,aes(x = fpr_100, y = tpr_100),color = "red")+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  annotate("text",x = .75, y = .25,
           label = paste("AUC of 60 = ",round(auc_60,2)),color = "blue")+
  annotate("text",x = .75, y = .15,label = paste("AUC of 100 = ",round(auc_100,2)),color = "red")+
  scale_x_continuous(name  = "fpr")+
  scale_y_continuous(name = "tpr")
```

