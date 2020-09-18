```R
options(stringsAsFactors = T)
sms_results <- read.csv('D:/download/sms_results.csv')
library(gmodels)
CrossTable(sms_results$actual_type,sms_results$predict_type)
library(caret)
#混淆矩阵
a <- confusionMatrix(sms_results$predict_type,
                sms_results$actual_type,
                positive = "spam")
#灵敏度： sensitivity=TP/(TP+FN)
sensitivity(sms_results$predict_type,
            sms_results$actual_type,
            positive = 'spam')
#特异度：specificity=TN/(TN+FP)
specificity(sms_results$predict_type,
            sms_results$actual_type,
            negative = 'ham')

library(ROCR)
data("ROCR.simple")
pred <- prediction(ROCR.simple$predictions, ROCR.simple$labels)
perf <- performance(pred,"tpr","fpr")
plot(perf,col="blue")
abline(a=0,b=1,lwd=2,lty=2)

###计算曲线下的AUC即面积
auc<-  performance(pred,"auc")
auc_area<-slot(auc,"y.values")[[1]] #auc_area<-auc@"y.values"[[1]]
auc_area<-round(auc_area,4)

#添加文本
text_auc<-paste("AUC=", auc_area,sep="")
text(0.8,0.3,text_auc)
```

