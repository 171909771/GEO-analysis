## 标准化矩阵
```
方法1
dge=edgeR::calcNormFactors(dge) ##先标准化
dge=cpm(dge,log = T) ##标准化矩阵，并且log2数值
方法2 直接看voom中的E值
```
