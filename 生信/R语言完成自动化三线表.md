## 三线表

```R
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
```

## 快速绘制临床论文基线特征表

```R
CreateTableOne(vars, strata, data, factorVars,
               includeNA = FALSE,
               test = TRUE,
               testApprox = chisq.test,
               argsApprox = list(correct = TRUE),
               testExact = fisher.test,
               argsExact = list(workspace = 2 * 10^5),
               testNormal = oneway.test,
               argsNormal = list(var.equal = TRUE),
               testNonNormal = kruskal.test,
               argsNonNormal = list(NULL),
               smd = TRUE, addOverall = FALSE)

# 函数参数解释说明
# vars  # 字符向量；指定哪些变量是基线特征表需要汇总的变量
# # 数据集中的因子视为分类变量，数字型变量视为连续变量
# # vars参数为空，则指定数据集中所有变量进行汇总
# 
# strata # 字符向量；指定分组汇总的变量，为空则进行单组汇总(也就是Overall列)
# data  # 变量来源的数据集名称，所有汇总变量都要在数据集里面
# factorVars # 字符向量；指定哪些变量为分类变量，指定的变量应是vars参数中的变量
# includeNA = FALSE # 逻辑词；为TRUE则将缺失值作为因子处理，仅对分类变量有效
# 
# test = TRUE # 逻辑词；默认为TRUE，当有2个或多个组时，自动进行组间比较
# testApprox = chisq.test # 默认卡方检验；当分类变量的单元格期望值较低(如<5）时，不建议使用
# argsApprox = list(correct = TRUE) # 进行chisq.test的连续校正。
# 
# testExact = fisher.test # 进行精确检验的函数，默认为fisher.test。
# argsExact = list(workspace = 2*10^5) # 指定fisher.test分配的内存空间
# 
# testNormal = oneway.test # 连续变量为正态分布进行的检验
# # 默认为oneway.test，两组时相当于t检验
# argsNormal = list(var.equal = TRUE) # 假设为等方差分析
# 
# testNonNormal = kruskal.test # 连续变量为非正态分布变量进行的检验，非参数检验
# # 默认为Kruskal-Wallis秩和检验
# # 当只有两个组时，与wilcox.test(Man-Whitney U检验)等效
# argsNonNormal = list(NULL) #传递给testNonNormal中指定的函数的参数的命名列表
# # 默认list(NULL)，它只是一个占位符
# 
# smd = TRUE # 默认为TRUE；当有两个以上的组时，则将自动计算组间比较的标准化均值差
# addOverall = FALSE # 仅在分组汇总中使用，将overall列添加到基线表中
# # smd和p值仅在分组汇总中使用。

library(tableone)
library(survival)  # 需要使用survival包的colon数据
data(colon)
CreateTableOne(data = colon)  # 汇总整个数据集特征
# 有两种方法可以将分类变量转化为因子：
# 一是先在数据集中将分类变量转化为因子，然后再使用tableone包进行汇总，
# 二是在tableone包中直接指定哪些变量属于因子(使用factorVars参数进行转换)，然后在进行汇总。
dput(names(colon)) # 输出colon数据集变量名称
#排除一些不需要比较的变量，如id、study之类的变量
myVars <- c("rx", "sex", "age", "obstruct", "perfor", 
            "adhere", "nodes", "status", "differ", 
            "extent", "surg", "node4", "time", "etype")
#指定基线表中哪些变量是分类变量
catVars <- c("rx", "sex","obstruct", "perfor", "adhere", "status",
             "differ","extent", "surg", "node4","etype")
tab2 <- CreateTableOne(vars = myVars,  
                       data = colon, 
                       factorVars = catVars)
# 通过vars参数指定哪些变量是基线表中需要汇总的变量
# 通过factorVars参数指定哪些变量是分类变量
# data参数指定变量的数据来源


#如果在表中要显示所有水平的数据，则输入：
print(tab2, showAllLevels = TRUE)
summary(tab2)

# 从前面的表中可以看到，连续变量都表示为均数＋标准差，这是认为连续变量都呈正态分布。
# 但是实际上有些数据呈非正态分布，需要用中位数（四分位数）表示。
#假设数据集中"time"和"nodes"两个连续变量呈非正态分布。
nonvar <- c("time","nodes") # 指定哪些变量是非正态分布变量
print(tab2, nonnormal = nonvar) # 输出基线表数据信息
#如果在print()函数中输入的是nonnormal = TRUE，则所有连续变量都按非正态分布进行分析。


tab3 <- CreateTableOne(vars = myVars, 
                       strata = "status", 
                       data = colon, 
                       factorVars = catVars)
# strata参数表示分层，指定需要分层的变量，这里我们指定status变量
# 通过vars参数指定哪些变量是基线表中需要汇总的变量。
# 通过factorVars参数指定哪些变量是分类变量
# data参数指定变量的数据来源


# CreateTableOne()函数默认的检验方法为：
# 分类变量使用卡方检验(chisq.test()，进行连续性校正)；
# 连续变量使用方差分析(oneway.test()，假设等方差)，两组间方差分析相当于t检验。

# 但是在基线表中，有些连续变量是非正态分布变量，
# 有些分类变量中单元格期望值较小，这些变量的统计方法不能使用默认的统计方法。
# 
# kruskal.test()函数可以用于呈非正态分布的连续变量，
# fisher.test()可以指定分类变量进行fisher精确检验。
# 在两组间比较时，kruskal.test()和wilcox.test()等效。

#基线表的test列会显示哪些变量是使用非默认检验方法来计算p值。

print(tab3, # 前面的tab3对象
      nonnormal = nonvar, # 指定哪些连续变量是非正态分布变量
      exact = "extent")  # 指定哪些变量需要使用fisher精确检验



tab4 <- CreateTableOne(vars = myVars, 
                       strata = "status", 
                       data = colon, 
                       factorVars = catVars, 
                       addOverall = TRUE) # 增加overall列
print(tab4, nonnormal = nonvar, exact = "extent")
# strata表示分层，指定需要分层的变量，这里我们指定status变量
# 通过vars参数指定哪些变量是基线表中需要汇总的变量。
# 通过factorVars参数指定哪些变量是分类变量
# data参数指定变量的数据来源


#输出基线特征表

# 简单粗暴的方法：就是复制粘贴，使用quote = TRUE显示引号，
# 使用noSpaces = TRUE删除用于在R控制台中对齐文本的空格，
# 然后直接复制基线表整个内容并将其粘贴到Excel电子表格即可。
print(tab4, # 前面的tab4对象
      nonnormal = nonvar, # 指定非正态分布变量
      exact = "extent", # 指定哪些变量需要使用fisher精确检验
      quote = TRUE,  # 显示引号
      noSpaces = TRUE) # 删除用于在R控制台中对齐文本的空格


#另一种方式：如果您不喜欢复制和粘贴，则可以通过以下方式自动导出。
tab4Mat <- print(tab4, nonnormal = nonvar, exact = "extent", 
                 quote = FALSE, # 不显示引号
                 noSpaces = TRUE, # 删除用于在R控制台中对齐文本的空格
                 printToggle = FALSE)

## 保存为 CSV 格式文件，并命名为 myTable。
write.csv(tab4Mat, file = "D:/linshiwenjian/r-daima/linshi/myTable.csv")


#仅输出分类变量
tab3$CatTable
#仅输出连续变量
print(tab3$ContTable, nonnormal = nonvar)
```

微调基线特征表输出格式。

```R
print(
  x,  # CreateTableOne函数创建的对象
  catDigits = 1, # 默认为1位，分类变量中百分比的小数位数
  contDigits = 2, # 默认为2位，连续变量的小数位数
  pDigits = 3, # 默认为3位，p值的小数位数，也可以设置标准化均值差异的小数位数
  quote = FALSE, # 默认为FALSE，是否在引号内显示内容；为TRUE，则将所有内容加上引号，以便复制到EXCEL
  missing = FALSE, # 默认为TRUE，是否显示缺失数据信息
  explain = TRUE, # 默认为TRUE，将“(%)”添加到分类变量名称中
  printToggle = TRUE, # 默认为TREU；为FALSE，则不创建任何输出
  test = TRUE, # 是否显示统计检验结果p值，默认为TRUE；为FALSE，则仅显示汇总结果
  smd = FALSE, # 是否显示标准化均值差异，默认FALSE
  noSpaces = FALSE, # 删除R控制平台中为对齐而产生的空格，复制到其他软件中时选TRUE。
  padColnames = FALSE, # 是否用空格填充列名以居中对齐，默认为FALSE，如果noSpaces = TRUE.则不进行
  varLabels = FALSE, # 是否用从labelled :: var_label（）函数获得的变量标签替换变量名。
  format = c("fp", "f", "p", "pf")[1], # 默认显示为“fp”频率（百分比），“f”为频率，“p”为百分比、“pf”为百分比（频率）。
  showAllLevels = FALSE, # 是否显示分类变量所有水平数据，默认为FALSE。
  cramVars = NULL, # 字符向量，指定二分类变量的两个水平在同一行显示
  dropEqual = FALSE, # 二分类变量中不显示第二水平名称；三分类及以上全都显示
  exact = NULL, # 字符向量，指定哪些分类变量进行精确检验，默认所有分类变量进行卡方检验
  nonnormal = NULL, # 字符向量，指定哪些连续变量进行非参数检验，默认所有连续变量为正态分布
  minMax = FALSE, # 逻辑词，默认为FALSE，非正态分布变量，是否使用[min,max]代替四分位数。
  ...
)
```

## t检验

计量资料的假设检验中，最简单、常用的方法就是 t 检验。

常见的t检验包括**单样本 t 检验，配对样本 t 检验和独立样本 t 检验**。

独立样本 t 检验一般要求数据服从**正态分布且方差齐性**。在进行 t 检验前一般先对资料进行方差齐性检验，若方差齐性，采用一般 t 检验；方差不齐，采用近似 t 检验。

配对样本 t 检验则要求每对数据差值的总体服从正态分布。

### 正态性检验

正态性检验的方法有两种，一是Q-Q图示法；二是正态 W 检验法。

#### Q-Q图示法

```R
qqnorm(x, ylim,
       main="NormalQ-QPlot",
       xlab="TheoreticalQuantiles",
       ylab="SampleQuantiles",
       plot.it=TRUE, datax=FALSE,...)
```

```R
options(digits = 3)  # 设定小数点后有效数字
x <- rnorm(20, mean=4, sd=4); x   # 生成正态分布数列
qqnorm(x); qqline(x)
```

若散点大致都在一条直线上，便可认为数据是服从正态分布的。

#### 正态W检验法

```R
H0:数据服从正态分布；
H1:数据不服从正态分布。
x <- rnorm(20, mean=4, sd=4); x  # 生成正态分布数列
shapiro.test(x)  # x是由数据构成的向量

#输出：
#Shapiro-Wilk normality test
#data:  x
#W = 1, p-value = 0.8
```

`p-value = 0.8 ＞ 0.05`，认为数据服从正态分布。

#### 方差齐性

方差齐性检验多采用 `Levene` 检验。

需要的包为 `car` 包，函数为`leveneTest()`函数

```R
library(car)
#leveneTest()函数接受数据框结构，一列是各分组的取值 y，另一列是分组 group。
#用法：leveneTest(y, group, center=median, ...)  # center可选 mean 和 median(默认)。

options(digits = 3)
x <- rnorm(20, mean=4, sd=4)# 生成数据
y <- rnorm(20, mean=5, sd=4)
d <- data.frame(x,y)  # 创建数据框
library(reshape2)
d1 <- melt(d, measure.vars = c("x","y"),  # 宽数据转长数据
           variable.name = "group",
           value.name = "groupvalue")

leveneTest(groupvalue ~ group, center = mean, data = d1) # 方差齐性检验
```

`P-value = 0.27 ＞ 0.05`，可认为等方差。

#### t 检验

```R
t.test(x,y = NULL,  # x，y是由数据构成向量（只提供x，做单个正太总体均值检验）
       alternative=c("two.sided","less","greater"),  # alternative是备择假设，括号内分别表示双尾(默认)、单尾（μ1＜μ2）和单尾（μ1＞μ2）检验
       mu=0,  # 表示原假设μ0
       paired=FALSE,  # 不配对，true为配对
       var.equal=FALSE,  # 默认方差不齐，TRUE表示方差齐性
       conf.level=0.95,...)  # 置信水平
```

##### 单样本 t 检验

```R
x <- rnorm(20, mean=5, sd=4); x

t.test(x, alternative="greater", mu=5) # 或 x-n  n为待检验均值

#输出：

#One Sample t-test
#data:  x
#t = -2, df = 19, p-value = 1
#alternative hypothesis: true mean is greater than 5
#95 percent confidence interval:
# 0.487   Inf
#sample estimates:
#mean of x
#     2.38 
```

##### 配对样本 t 检验

```R
x <- rnorm(20, mean=5, sd=4); x
y <- rnorm(20, mean=4, sd=4); y

t.test(x-y, alternative="two.sided")

# 输出：
# 
# One Sample t-test
# data:  x - y
# t = 1, df = 19, p-value = 0.2
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   -0.621  3.585
# sample estimates:
#   mean of x
# 1.48 
```

##### 独立样本 t 检验

```R
# 生成两组数据
x <- rnorm(20, mean=5, sd=4); x
y <- rnorm(20, mean=4, sd=4); y

t.test(x,y, var.equal=TRUE, alternative="two.sided") # 独立样本t检验

# 输出：
# 
# Two Sample t-test
# data:  x and y
# t = 1, df = 38, p-value = 0.2
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.562  3.526
# sample estimates:
#   mean of x mean of y
# 4.21      2.73
```

若方差不齐，不能满足 t 检验，采用 `Wilcoxon` 秩和检验，也叫 `Mann-Whitney` 检验。

```R
wilcox.test(y ~ x, data)  # y是数值型变量，x为二分类变量，data为矩阵或数据框 或
wilcox.test(y1, y2) # y1、y2为各组的数值型向量

# wilcox.test(x, y = NULL,
#             alternative = c("two.sided", "less", "greater"),
#             mu = 0,
#             paired = FALSE,
#             exact = NULL,  # 逻辑词，是否计算精确p值
#             correct = TRUE,
#             conf.int = FALSE, conf.level = 0.95, ...)
```

