## 标准化矩阵
```
dge=edgeR::calcNormFactors(dge) ##先标准化
dge=cpm(dge,log = T) ##标准化矩阵，并且log2数值
```
