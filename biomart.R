# <差异基因分析>
# 1.判断是否有BiocManager包，若不存在则安装
#options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))) #设置清华镜像，加速下载

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#if (!requireNamespace("DESeq2", quietly = TRUE))
#  BiocManager::install('DESeq2')  #通过BiocManager安装DESeq2 
library(DESeq2) #加载library

setwd("C:\\Users\\lzhang\\Desktop\\TGEV") #设置工作目录，所有输出文件保存于此

#输入数据要求
# DEseq2要求输入数据是由整数组成的矩阵
# DESeq2要求矩阵是没有标准化的

##2.读入所有基因原始readscount表达矩阵，行为基因，列为样品

A <- read.table("C:\\Users\\lzhang\\Desktop\\TGEV\\rsem.merged.gene_counts.tsv", header = T, row.names = 1)
B <- as.matrix(A) #转换成矩阵格式，保证都是数值

#View(B)
## 3.实验分组
# 样品信息矩阵即上述代码中的colData，它的类型是一个dataframe（数据框），
# 第一列是样品名称，第二列是样品的处理情况（对照还是处理等），即condition
coldata <- read.table("sample_info.txt",header = T,row.names = 1)
coldata <- coldata[, c("condition", "type")]
#View(coldata)	#查看分组信息

## 4.制作dds对象，构建差异基因分析所需的数据格式
dds <- DESeqDataSetFromMatrix(countData = B, colData = coldata, design = ~ condition);

# countData = B，readscount矩阵
# colData = coldata,分组信息，根据这个才能在2组之间比较
# design = ~ condition，公式，表示按照condition进行分析

## 5.差异分析结果
dds <- DESeq(dds)	#正式进行差异分析

## 6.提取结果，在treated和untreated组进行比较
res <- results(dds, contrast = c("condition", "WT", "PI")) 
# results从DESeq分析中提取出一个结果表，从而给出样品的基本均值，log2倍变化，标准误差，测试统计量，p值和校整后的p值； 
sum(res$padj < 0.05, na.rm = TRUE)	#统计padj小于0.05显著差异的基因

##8.过滤上调、下调基因
filter_up <- subset(res, pvalue < 0.05 & log2FoldChange > 1) #过滤上调基因
filter_down <- subset(res, pvalue < 0.05 & log2FoldChange < -1) #过滤下调基因
filter_diff  <- subset(res, padj < 0.05)	#统计padj小于0.05显著差异的基因

print(paste('差异上调基因数量: ', nrow(filter_up)))  #打印上调基因数量
print(paste('差异下调基因数量: ', nrow(filter_down)))  #打印下调基因数量

##9.保存到文件
write.table(filter_diff, file = "./differential_gene.txt", sep = "\t") #log2FoldChange + pvalue + padj
write.table(filter_up, file="./filter_up_gene.txt", quote = F, sep = "\t")  
write.table(filter_down, file="./filter_down_gene.txt", quote = F, sep = "\t")


#------------------------------------------------------
library(EnhancedVolcano)
#devtools::install_github('kevinblighe/EnhancedVolcano')
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-20,20),
                title ='resistant versus control',
                pCutoff = 10e-17,
                FCcutoff = 1 ,
                colAlpha = 1,
                col=c(' black','blue',' green','red1'),
                
)

#----------------------------------------------------------

library("biomaRt")
mart <- useDataset("sscrofa_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id <- row.names(filter_diff)
head(my_ensembl_gene_id)
pig_symbols <- getBM(attributes = c('ensembl_gene_id','external_gene_name','description'),filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(pig_symbols)
ensembl_gene_id <- rownames(filter_diff)
filter_diff <- cbind(ensembl_gene_id,filter_diff)
colnames(filter_diff)[1]<-c("ensembl_gene_id")
diff_name <-merge(filter_diff, pig_symbols, by="ensembl_gene_id")
diff_name
write.table(diff_name, file="./all_diff_genename.txt", quote = F, sep = "\t")  

