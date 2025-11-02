# Programming-Project
Programming Project

我打算先用一个小样本（比如说某个组织，甲状腺...）跑通这个pipeline，
然后通过将所有的组织样本跑通，合并在一起进行后续分析。

RNA seq数据：
https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
目前先采用了Gene read counts 中thyroid的部分

annotation（用于TPM正则化）：
https://www.gencodegenes.org/human/release_44.html
Comprehensive gene annotation	CHR（基因长度注释）

目前进度:
阶段1：数据清洗
a. 读取GTF: gtf_parser.go 负责解析 gencode.v44.annotation.gtf.gz 文件，计算出每个基因的准确长度
b. 加工数据: gct_processor.go 利用 gtf_parser.go 拿到的基因长度，对gene_reads_v10_thyroid.gct.gz执行一个“两遍流式处理”： 
  i. 第一遍: 计算TPM标准化所需要的分母（perSampleRPKSum）
  ii. 第二遍: 计算 TPM -> Log2(TPM+1) 转换
  iii. 两次过滤： “低表达过滤” -> “低变异过滤”
c. 产出: main.go 把这个做好的数据矩阵保存为 clean_thyroid_matrix.csv
