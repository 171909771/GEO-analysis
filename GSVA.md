library(GSVA)
library(GSEABase)
library(limma)

# 读取数据
## txt.gz 用read.table
```r
dat1=read.table("GSE22255_series_matrix.txt.gz",comment.char = "!",header = TRUE,row.names  =1)
```
# 转换探针名字
```r
library("hgu133plus2.db")
ids=toTable(hgu133plus2SYMBOL)
```r
##2 download GPL soft infromation from GEO and need some shell command to reconstructed data
## keep repeated symbl by add number post the names
```r
if(F){
  ids$symbol=make.unique(ids$symbol)   ###使重复基因添加序号
  dat1=as.data.frame(dat1)
  dat1=dat1[ids$probe_id,]
  dat1$symbl=ids$symbol
  rownames(dat1)=dat1$symbl
  dat1$symbl=NULL
}
## unique the rownames
dat1=as.data.frame(dat1)
dat1$md=apply(dat1,1,sd)
dat1=dat1[ids$probe_id,]
dat1$symbl=ids$symbol
dat1=dat1[with(dat1,order(symbl,md,decreasing = T)),]
dat1=dat1[!duplicated(dat1$symbl),]
rownames(dat1)=dat1$symbl
dat1$symbl=NULL
dat1$md=NULL

## 转换ENTREZID
library(clusterProfiler)
gene <- bitr(rownames(dat1),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = "org.Hs.eg.db")
dat1$genename=rownames(dat1)
dat2=merge(dat1,gene,by.x="genename",by.y="SYMBOL")
rownames(dat2)=dat2$ENTREZID
dat2$genename=NULL
dat2$ENTREZID=NULL
```
# GSVA
## 转换表达矩阵到通路矩阵
```r
keggSet <- getGmt("c2.cp.kegg.v7.5.1.entrez.gmt")
keggEs <- gsva(expr=as.matrix(dat2), gset.idx.list=keggSet, kcdf="Poisson", parallel.sz=4)
head(keggEs, n=3)
```
## 建立对比关系
```r
grouP <- c(rep("A", 20), rep("B", 20)) %>% as.factor()
desigN <- model.matrix(~ grouP + 0)
rownames(desigN) <- c(paste("A",1:20,sep = ""),paste("B",1:20,sep = ""))
desigN
### 用desigN的列名，设定B组比A组
comparE <- makeContrasts(grouPB - grouPA, levels=desigN)
```
## 差异分析
```r
fiT <- lmFit(keggEs, desigN)
fiT2 <- contrasts.fit(fiT, comparE)
fiT3 <- eBayes(fiT2)
keggDiff <- topTable(fiT3, coef=1, number=200)
head(keggDiff, n=3)
```
