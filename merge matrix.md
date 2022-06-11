## https://cloud.tencent.com/developer/article/1587442

  ## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48213
  #wget -c ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48213/suppl/GSE48213_RAW.tar
  #tar -xf GSE48213_RAW.tar
  #gzip -d *.gz
  ## 首先在GSE48213_RAW目录里面生成tmp.txt文件，使用shell脚本
  # awk '{print FILENAME"\t"$0}' GSM*.txt |grep -v EnsEMBL_Gene_ID >tmp.txt
  #  其实也可以直接使用R来读取GSE48213_RAW.tar里面的gz文件，这里就不演示了
  # 可以参考：https://mp.weixin.qq.com/s/OLc9QmfN0YcT548VAYgOPA 里面的教程
  ## 然后把tmp.txt导入R语言里面用reshape2处理即可
  # 这个 tmp.txt 文件应该是100M左右大小哦。
  
```r
  a=read.table('GSE48213_RAW/tmp.txt',sep = '\t',stringsAsFactors = F)
  library(reshape2)
  fpkm <- dcast(a,formula = V2~V1) 
```
## 合并多个csv
```r
dat=list.files("./")
data=data.frame(name=read.table(path, header = TRUE)[1])
### 方法1，更快，但是要确定合并的名字相同
for (i in dat){
  path=i
  data <- cbind(data, read.table(path, header = TRUE)[2])
}
### 方法2：更慢，可以按照相同的名字合并
for (i in dat){
  path=i
  data <- merge(data, read.table(path, header = TRUE),by = "Geneid")
}
```
### 样式
![image](https://user-images.githubusercontent.com/41554601/173173096-6e480272-c897-41b5-8fbf-42b7a4805548.png)
