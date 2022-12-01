# edgeR、limma、DESeq2三种差异表达包比较（RNA-seq数据）

- https://blog.csdn.net/weixin_43700050/article/details/98085127
## limma包的详细应用
- https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow_CHN.html

# limma处理基因芯片数据，主要是用RMA归一化
- https://liripo.netlify.app/post/%E5%9F%BA%E5%9B%A0%E8%8A%AF%E7%89%87%E5%B7%AE%E5%BC%82%E6%80%A7%E5%88%86%E6%9E%90/#fn:2

# GEO FPKM 转 TPM后limma分析 + 富集网络图
- https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247513130&idx=4&sn=622366cd28226a70a8d3961e5ffad0a6&chksm=9b4bf491ac3c7d87f2a21817b45419e0f084c340692b56e0edaf55b78c83422811a7656d762b&mpshare=1&scene=1&srcid=0423ZyMUSNQnrsnJcvgeVSJE&sharer_sharetime=1650698239528&sharer_shareid=0cdcf284bc426c9af2f9246b33b7ac7f&exportkey=Aaqj1anPf4%2FgrfW%2BN96VYkA%3D&acctmode=0&pass_ticket=2QHlunDxvd1colkYwJ5ehwY2emfc3KhhlU1BsrZnZikJRIPDXZoPP%2BEii5q0r9GB&wx_header=0#rd

# 当每个条件的样本量大于8个以后，别犹豫，直接上Wilcoxon秩和检验。
- https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247512042&idx=1&sn=acf55e1a4090661b75dca179213ecfd5&chksm=9b4bef51ac3c6647e387e8c9ab057430ea2c50487babf1f0eef7a62cdb7194d4f6a5077699db&mpshare=1&scene=24&srcid=0331JBdCxi66Ja7KQMqpBruA&sharer_sharetime=1648731156544&sharer_shareid=0cdcf284bc426c9af2f9246b33b7ac7f&exportkey=AdTeJnb%2FSpNK0OeKcQZULIY%3D&acctmode=0&pass_ticket=DAcOzZoucYKsK8h7SbdTKTlVTWXufLfNse3D5orsWIXSt0%2BJDXeZNhQX3EIDP5sW&wx_header=0#rd

# miRNA 和 mRNA的关联分析
## step1 先计算miRNA的差异基因（可以用edgeR、limma、DESeq2来分析）
## step2 在把数据放到miRWalk 3.0 和 miRTarBase databases中分析，出关联图
### miRWalk 3.0可能需要自己再分析一次

# rpkm转TPM
```r
library("barzinePhdR")
exp1=rpkm2tpm(exp)
```

# edgeR归一化TMM方法详解
- https://www.jianshu.com/p/935a5466b8c2

# 讲解DESeq2 的原理，其中的dispersion非常精彩
- https://wap.sciencenet.cn/blog-3431904-1308379.html?mobile=1
##### 辅助讲解负二项分布中dispersion的由来
- https://zhuanlan.zhihu.com/p/111632687
#### 需要看下里面的PCA 和聚类的连接
- https://cloud.tencent.com/developer/article/1949605
