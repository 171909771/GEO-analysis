load("dat1.Rdata")

library(WGCNA)

# 前处理数据
## 读取矩阵，已经FPKM或则其他类似的标准化的矩阵
```r
datExpr <- t(keggEs)  #转置
```
## 检查样本的聚类情况，有没有异常组
```r
datExpr0=datExpr
sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "sampleClustering.pdf", width = 12, height = 9)
par(cex = 1)
par(mar = c(0,3,4,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 90, col = "red") #先画一条辅助线
![1](https://user-images.githubusercontent.com/41554601/170813414-4189860d-d184-4248-a2eb-90e4b4ba7246.png)
## 聚类后裁枝
clust = cutreeStatic(sampleTree, cutHeight = 90, minSize = 10)
table(clust) # 0代表切除的，1代表保留的
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
```

# 正式开始WGCNA
## 计算软阈值
### 可以把R方设置到0.8
```r
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
best_beta=sft$powerEstimate
```

## 计算相关性，建立modules
```r
cor <- WGCNA::cor
net = blockwiseModules(
  datExpr,
  power = 6,
  maxBlockSize = 5000, #default is 5000
  TOMType = "unsigned", minModuleSize = 30,   
  reassignThreshold = 0, mergeCutHeight = 0.25,  # mergeCutHeight: 合并模块的阈值，越大模块越少；越小模块越多，冗余度越大；
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = F, 
  verbose = 3
)
table(net$colors)

### 每个模块上色
mergedColors = labels2colors(net$colors)
table(mergedColors)
### 根据模块颜色绘图
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
```

## 最重要的模块与表型的相关图
```r
- https://www.biostars.org/p/414257/    ## 多级分类性状的处理  
## 这一步主要是针对于连续变量，如果是分类变量，变换成0，1代表，如果是多分类，变换成0，1，3......，也可以按照天数变换0.3.5.6、、、、
示例数据
[meta.csv](https://github.com/171909771/GEO-analysis/files/8790608/meta.csv)
```r
#### 方法一：计算单个多分类变量
datTraits=model.matrix(~0+ datTraits$subtype)
colnames(datTraits)=levels(as.factor(datTraits$subtype))
#### 方法二：计算多个多分类变量
datTraits=read.csv("meta1.csv")
rownames(datTraits)=datTraits$ID 
datTraits$ID=NULL
```
### 开始计算模块与表型之间的相关性
```r
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#### 用color labels重新计算MEs（Module Eigengenes:模块的第一主成分）
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
#### 排序
MEs = orderMEs(MEs0)
####（这是重点）计算ME和表型相关性
moduleTraitCor = cor(MEs, datTraits, use = "p") 
```
### 出图
```r
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
```


## 确定目标模块中的基因在性状和模块中都具有意义
```r
  # 先把模块的颜色名字去掉前面的“ME”
  modNames = substring(names(MEs), 3)
  # 算出每个模块跟基因的皮尔森相关系数矩
  ## MEs是每个模块在每个样本里面的相关系数矩
  ## datExpr是每个基因在每个样本的表达量
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  names(geneModuleMembership) = paste("MM", modNames, sep="");


  # 计算基因在性状上的相关系数
  Luminal = as.data.frame(datTraits[,1])
  names(Luminal) = "Luminal"
  geneTraitSignificance = as.data.frame(cor(datExpr, Luminal, use = "p"));
  names(geneTraitSignificance) = paste("GS.", names(Luminal), sep="");

  # 整合上面的数据出图
  module = "white"
  column = match(module, modNames);
  moduleGenes = mergedColors==module;
  png("step6-Module_membership-gene_significance.png",width = 800,height = 600)
  #sizeGrWindow(7, 7);
  par(mfrow = c(1,1))
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for Luminal",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "red")
```
  


