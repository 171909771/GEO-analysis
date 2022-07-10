
- https://zhuanlan.zhihu.com/p/115305225
- https://cloud.tencent.com/developer/article/1643009

## 数据前处理
```
library(GEOquery)
library(Biobase)
library(limma)
library(AnnoProbe)
### 下载raw数据，以GSE110141为例，下载并且解压，示例数据文件夹为“sampFile”
dat = list.files("sampFile/", pattern = "txt")
### agilent 特有的读取方法
dat=read.maimages(files= dat,path = "./sampFile",green.only = TRUE,source="agilent")
### 背景校正和标准化处理
RG = backgroundCorrect(dat, "normexp", normexp.method = "rma", offset=50) 
E = normalizeBetweenArrays(RG, method="quantile")
### 构建表达矩阵
ep=as.data.frame(E)
ep=ep[,6:length(colnames(ep))]
### 看标准化后的boxplot水平
#### rainbow是R的一个函数，用于产生彩虹色
par(cex = 0.9)
cols <- rainbow(length(ep) * 1.2)
boxplot(ep, col = cols, xlab = "Sample", ylab = "Log intensity")
### 加入probeid
ep$probe_id <- E[["genes"]]$ProbeName
### 相同的探针取平均值
ep <- aggregate(x = ep[,1:12],
                by = list(ep$probe_id),
                FUN = mean)
```

## 基因注释
```
### 读取注释信息
gpl=idmap(gpl = "GPL10787", type = "pipe", mirror = "tercent")
if(F){ ### 多个芯片平台
gpl1=idmap(gpl = "GPL10787", type = "pipe", mirror = "tercent")
gpl2=idmap(gpl = "GPL21163", type = "pipe", mirror = "tercent")
gpl=unique(rbind(gpl1,gpl2))
}
### 合并data和注释信息
dat1=merge(ep,gpl,by.x ="Group.1",by.y ="probe_id")
### 相同的基因取平均值
expr <- aggregate(x = dat1[,2:13],
                      by = list(dat1$symbol),
                      FUN = mean)
### 生成最终的表达矩阵
rownames(expr)=expr$Group.1
expr$Group.1=NULL
```
