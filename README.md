# test1
|                                                              |
| ------------------------------------------------------------ |
| 从TCGA下载Tissue Slide数据<br/>/home/zhoulin/gdc-client download -m coad.txt(manifest文件) |
| 生成切片<br/>conda activate os<br/>nohup bash coad.sh &<br/>nohup bash test.sh >> test.txt 2>&1 & |
| anconda环境<br/>os                       /home/zhoulin/.conda/envs/os: openslide<br/>pytor                 *  /home/zhoulin/.conda/envs/pytor：torch<br/>tfg                      /home/zhoulin/.conda/envs/tfg |
| editplus:注释--》ctrl+d  取消注释--》ctrl+b                  |
| `trainCoxMlp` -训练Cox-nnet模型的主要功能 `CVLoglikelihood` -计算交叉验证的对数可能性（模型性能指标） `CIndex` -计算C-Index（模型性能指标） `L2CVSearch`-用于优化正则化参数的辅助函数；使用爬山算法搜索lambda `L2CVProfile`-用于优化正则化参数的辅助函数；跨一系列值配置lambda `evalNewData`-评估新数据或测试数据；输出线性预测因子（即对数风险比） `varImportance` -通过辍学程序确定变量的重要性 `saveModel` -将模型保存到二进制文件 `loadModel` -从文件加载模型 |

It's a test repository

./bget.exe doi 10.1186/s12935-020-01351-3 -t 5 --suppl

![Error](https://github.com/bigone1/test/blob/master/Screenshots/1.png)

| Git 本地上传指令                                   |
| -------------------------------------------------- |
| git add .                                          |
| git commit -m "tj"                                 |
| git pull origin master --allow-unrelated-histories |
| git push origin master                             |

=VLOOKUP(A2,Sheet1!A:M,10,0)

插件下载：Crx4Chrome

| 网站   | 账户                    | code       |
| ------ | ----------------------- | ---------- |
| Oracle | zl0508@mail.ustc.edu.cn | Zl159357+- |
| 58     | 18755158902             | zl159357   |
|        |                         |            |

```R
withr::with_makevars(c(PKG_LIBS = "-liconv"), install.packages("haven"), assignment = "+=")
withr::with_makevars(c(PKG_LIBS = "-liconv"), install.packages("readxl"), assignment = "+=")
install.packages("tidyverse")

for(i in 1:(nrow(g)-1)){
  print(i)
  for(j in (i+1):nrow(g)){
    pair <- ifelse(g[i,] > g[j,],1,0)
    pairRatio <- sum(pair)/sampleNum
    if((pairRatio > 0.2) & (pairRatio < 0.8)){
      rownames(pair) <- paste0(rownames(g)[i],"|",rownames(g)[j])
      f <- rbind(f,pair)
    }
  }
}
```

