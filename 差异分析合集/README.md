## 归一化总结
- https://blog.csdn.net/weixin_46585008/article/details/109478030
![image](https://user-images.githubusercontent.com/41554601/200205289-0384f0d8-3f65-454d-a4c8-73703e24454e.png)

## 标准化矩阵
```
方法1
dge=edgeR::calcNormFactors(dge) ##先标准化
dge=cpm(dge,log = T) ##标准化矩阵，并且log2数值
方法2 直接看voom中的E值
```
