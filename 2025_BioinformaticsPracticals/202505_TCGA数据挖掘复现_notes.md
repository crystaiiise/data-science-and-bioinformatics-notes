**b站_10小时TCGA数据挖掘及文章复现**

Course link: <https://www.bilibili.com/video/BV1nc411g7G3>

// 本来是在R markdown里面跑来着，但反正最后也没knit几个代码块，不如直接拷贝到md里面方便。

---

- [Day 1: Basic R Recaps](#day-1-basic-r-recaps)
    - [Recap: data.frame](#recap-dataframe)
  - [Recap: processing data](#recap-processing-data)
- [Day 2: 访问并整理TCGA数据](#day-2-访问并整理tcga数据)
  - [访问](#访问)
  - [数据提取：counts及TPM](#数据提取counts及tpm)
  - [根据样本Barcode整理数据](#根据样本barcode整理数据)
- [Day 3: 免疫得分计算与生存分析](#day-3-免疫得分计算与生存分析)
  - [几种得分计算](#几种得分计算)
  - [KM生存分析](#km生存分析)
- [Day 4: 临床数据整理及绘图](#day-4-临床数据整理及绘图)
- [Day 5: 差异分析](#day-5-差异分析)
    - [理解、定义并筛选差异基因](#理解定义并筛选差异基因)
    - [热图绘制及表达谱样本整理](#热图绘制及表达谱样本整理)
  - [两次差异分析结果取交集](#两次差异分析结果取交集)
  - [用交集差异基因做GO和KEGG富集分析](#用交集差异基因做go和kegg富集分析)
    - [富集分析后绘制网络图](#富集分析后绘制网络图)
- [Day 6: 蛋白质互作及COX回归分析](#day-6-蛋白质互作及cox回归分析)
  - [PPI](#ppi)
    - [STRING与Cytoscape](#string与cytoscape)
  - [单因素COX回归分析](#单因素cox回归分析)
  - [森林图绘制及读图](#森林图绘制及读图)
    - [解决图像显示/导出不完整问题！！](#解决图像显示导出不完整问题)
  - [取PPI与Cox分别筛选后的基因交集](#取ppi与cox分别筛选后的基因交集)
- [Day 7-8: 用筛选到的单基因BTK做各种分析](#day-7-8-用筛选到的单基因btk做各种分析)
  - [正常/肿瘤样本中表达量比较-柱状图](#正常肿瘤样本中表达量比较-柱状图)
  - [样本配对原理及配对图](#样本配对原理及配对图)
  - [用单基因表达量分高低组做生存分析](#用单基因表达量分高低组做生存分析)
  - [比较不同临床分期的样本中BTK表达量](#比较不同临床分期的样本中btk表达量)
  - [BTK表达量高低组的差异分析及富集分析](#btk表达量高低组的差异分析及富集分析)
    - [用到本地基因集gmt注释文件的GSEA](#用到本地基因集gmt注释文件的gsea)
    - [富集分析绘图并读图](#富集分析绘图并读图)
- [Day 9: CIBERSORT免疫浸润分析](#day-9-cibersort免疫浸润分析)
  - [彩虹图绘制及读图](#彩虹图绘制及读图)
  - [相关性热图](#相关性热图)
  - [BTK表达量高低组的免疫细胞分组比较图](#btk表达量高低组的免疫细胞分组比较图)
  - [BTK单基因表达量与免疫细胞占比 相关性散点图](#btk单基因表达量与免疫细胞占比-相关性散点图)
    - [一定区分上述的差异分析与相关性分析](#一定区分上述的差异分析与相关性分析)
- [总结：文章讲了什么故事...?](#总结文章讲了什么故事)
- [完结撒花！！！！](#完结撒花)


## Day 1: Basic R Recaps

新建 "New Project"，这样子每次通过选定文件夹里的.Rproj文件打开喵。

（从Day9回来，想找个地方着重记一下：真正到画图等等的步骤，数据整理的成果要让**想研究的变量作为输入数据中的列**。）

#### Recap: data.frame 

```R
df <- read.table("x.txt", sep = "\t", row.names = 1, check.names = F, header = T)
write.csv(df, file = "day1.csv")

# Transposing a df
t_df <- t(df)
class(t_df) # so we would have to convert it back to a df
df <- as.data.frame(t_df)
# Equivalent to (using pipeline):
library(tidyverse)
df <- df %>% t() %>% as.data.frame()
df
```

### Recap: processing data

这个其实是Day2涉及的但是就放在这边吧qvq

```R
library(tidyverse)
# Remove duplicates
a <- c(1,1,1,1,1,3,4,5)
duplicated(a)
a[!duplicated(a)] # Subsetting by logical vector


# Tibble recap 
# Note that tribble() is just a function that refers to "transposed tibble", allowing us to create tibbles in a more readable, **row-wise** manner!)
t1 <- tribble(
  ~"col1",~"col2", # Set the variables of the (transposed) tibble; fixed!
  "a",1L,
  "b",2L
)
# This is just equivalent to the ordinary tibble() by:
t2 <- tibble("col1" = letters[1:2], "col2" = 1:2) # Column-wise
t1
t1 == t2
identical(t1,t2) # (Super strict, if we didn't specify the numeric values in t1 to be integers (as returned by t2's `1:2`), the output of identical() will still be FALSE)
t3 <- tibble("col1" = letters[1:4], "col3" = 6:9)
inner_join(t1,t3,by = "col1")
right_join(t1,t3,by = "col1")
```

---

## Day 2: 访问并整理TCGA数据

### 访问

```R
dir.create("TCGA-LUAD/TCGAdata", recursive = T)
setwd("TCGA-LUAD/TCGAdata")
library(BiocManager)
library(TCGAbiolinks)
```

TCGA数据库中癌症缩写：<https://www.jianshu.com/p/3c0f74e85825>。比如我们这里的LUAD就是Lung adenocarcinoma。

`GDCdownload()`中各参数：<https://zhuanlan.zhihu.com/p/80426727>。


```R
cancer_type <- "TCGA-LUAD"
expquery <- GDCquery(project = cancer_type,
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "STAR - Counts")
                     # "STAR" workflow: counts/FPKM/TPM data
GDCdownload(query = expquery, directory = "TCGA-LUAD/TCGAdata/GDCdata") 

# GDCprepare(): reads data into an R object

expquery2 <- GDCprepare(query = expquery, directory = "TCGA-LUAD/TCGAdata/GDCdata", summarizedExperiment = T) 
save(expquery2,file = "luad.gdc_2022.rda")
```

实在太慢了下载不下来，直接导入（.rda的话也可直接双击文件名）从课程文件夹下载的代码：
```R
setwd("TCGA-LUAD/TCGAdata")
load("luad.gdc_2022.rda")
load("gene_annotation_2022.rda")
table(gene_annotation_2022$type) # count the occurrence of each value in 'type' column
```

我们之前不是有这个嘛：（虽然没跑）
```R
expquery2 <- GDCprepare(query = expquery, directory = "TCGA-LUAD/TCGAdata/GDCdata", summarizedExperiment = T) 
save(expquery2,file = "luad.gdc_2022.rda")
```

`luad.gdc_2022.rda`导入进来后环境里得到`expquery2`。打开看看喵！——发现这个expquery2是**SummarizedExperiment 对象**（参数中有标明）。取measurement matrices可不就是用`@`取slot `assays`里面的嘛！

别忘了这些表达矩阵都是行作为features/基因，列为样本。


> Recalling JHU Course 5... 

（没想到十年前的课里面讲的概念还是这么有帮助的说  实战里面才开始串起来了OvO！） "The measurement data is accessed by `assay` and `assays`. A SummarizedExperiment can contain multiple measurement matrices (all of the same dimension!). "

```R
expquery2

# Output: 
# class: RangedSummarizedExperiment 
# dim: 60660 600 
# metadata(1): data_release
# assays(6): unstranded stranded_first ...
#   fpkm_unstrand fpkm_uq_unstrand
# rownames(60660): ENSG00000000003.15
#   ENSG00000000005.6 ... ENSG00000288674.1
#   ENSG00000288675.1
# rowData names(10): source type ...
#   hgnc_id havana_gene
# colnames(600):
#   TCGA-44-8120-01A-11R-2241-07
# ...
# ...

```



### 数据提取：counts及TPM

我们之前GDCquery的时候不是`workflow.type = "STAR - Counts"`嘛，除了counts以外实际下载的是STAR workflow生成的标准化数据，包含多种标准化格式（FPKM/TPM），所以能直接提取TPM不用从FPKM转换。

注：差异分析只能用counts做（这里称"unstranded"），毕竟是要raw data，而TPM是已经标准化过后的了；除了差异分析用counts之外其他各种下游分析都用其他measurement matrices（eg. TPM）做。

那么我们就来提取counts（记得S4类中的slots是用`@`索引）。提取完了之后利用基因注释那个数据给counts数据去重，可以发现counts从原本的六万多剩下了五万多行～～

```R
# Extract the counts data into a matrix and then name the columns and rows. (The following three lines are powerful in extracting data from any queried TCGA datasets; we do the same thing below when extracting TPM data)

# Recall that data inside an S4 class (just like our RangedSummarizedExperiment class here) is organized into slots. We access slots by using the @ operator. 
counts <- expquery2@assays@data@listData[["unstranded"]]
colnames(counts) <- expquery2@colData@rownames
rownames(counts) <- expquery2@rowRanges@ranges@NAMES

# Filtering the counts data to remove duplicated records (using ENSEMBL IDs from the gene annotation data)
counts <- counts %>% 
  as.data.frame() %>% 
  rownames_to_column("ENSEMBL") %>% 
  inner_join(gene_annotation_2022,"ENSEMBL") %>% 
  .[!duplicated(.$symbol),] # Note that you put a `.` to represent the referenced object in a pipeline!!!!
  # column_to_rownames(var = "symbol") note that this function requires the df to have no row names (even not the numeric indices) so we have to set it to NULL.
rownames(counts) <- NULL
counts <- counts %>% column_to_rownames(var = "symbol")

```


### 根据样本Barcode整理数据

记得上面我们用table()看过每条基因的type字段都是什么值/出现了多少次。现在我们只要mRNA，即只保留基因类型为"protein_coding"的行。之后对样本（列）进行去重，再提出来特定conditions下的样本。

列名称Barcode，有特定命名方法：[一文讲清TCGA数据库中样本编码信息](https://zhuanlan.zhihu.com/p/564801425)。注意看Sample和Vial，比如Sample中01为primary solid tumor, 11为solid tissue normal。我们在这里做差异分析看的就是01A（癌症组织）与11A（正常组织）。
```R
counts <- counts[counts$type == "protein_coding",]

# After removal of duplicates and gene-type filtering we no longer need the first column (ENSEMBL ID) and the last column (type). Negative indexing!
counts <- counts[, -c(1, ncol(counts))] 

colnames(counts) <- substring(colnames(counts),1,16)
counts <- counts[,!duplicated(colnames(counts))]

# Viewing the conditions of the samples (substrings for conditions (eg. '01A') now occur on the 14-16th characters of the entire column name)
table(substring(colnames(counts),14,16))
# Extracting
counts01A <- counts[,substring(colnames(counts),14,16) == "01A"]
counts11A <- counts[,substring(colnames(counts),14,16) == "11A"]
```

对于TPM，过程跟counts基本一模一样。TPM为标准化之后的，所以变成了连续性的，点开看看里面有小数（不像counts全部是整数，毕竟那个是reads计数嘛），有的太大了所以取个对数。
```R
tpms <- expquery2@assays@data@listData[["tpm_unstrand"]]
colnames(tpms) <- expquery2@colData@rownames
rownames(tpms) <- expquery2@rowRanges@ranges@NAMES
tpms <- tpms %>% 
  as.data.frame() %>% 
  rownames_to_column("ENSEMBL") %>% 
  inner_join(gene_annotation_2022,"ENSEMBL") %>% 
  .[!duplicated(.$symbol),]
rownames(tpms) <- NULL
tpms <- tpms %>% column_to_rownames("symbol") 
tpms <- tpms[tpms$type == "protein_coding",]
tpms <- tpms[,-c(1,ncol(tpms))]
colnames(tpms) <- substring(colnames(tpms),1,16)
tpms <- tpms[,!duplicated(colnames(tpms))]
tpms01A <- tpms[,substring(colnames(tpms),14,16) == c("01A")]
tpms11A <- tpms[,substring(colnames(tpms),14,16) == c("11A")]

range(tpms) # Massive; we apply the log transformation
tpms01A_log2 <- log2(tpms01A+1)
range(tpms01A_log2)
tpms11A_log2 <- log2(tpms11A+1)
range(tpms11A_log2)
```

判断counts和TPM的行列名是否一致，并保存数据喵。

```R
identical(rownames(counts01A),rownames(counts11A))
identical(rownames(tpms01A),rownames(tpms11A))
identical(rownames(counts01A),rownames(tpms01A))
identical(colnames(counts01A),colnames(tpms01A))

write.table(counts01A,"TCGA-LUAD/TCGAdata/counts01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(counts11A,"TCGA-LUAD/TCGAdata/counts11A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms01A,"TCGA-LUAD/TCGAdata/tpms01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms11A,"TCGA-LUAD/TCGAdata/tpms11A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms01A_log2,"TCGA-LUAD/TCGAdata/tpms01A_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms11A_log2,"TCGA-LUAD/TCGAdata/tpms11A_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
```
确认完了行名一致，我们可以用`cbind()`按列把01A和11A两种条件的样本拼起来。
```R
counts <- cbind(counts01A,counts11A)
table(substring(colnames(counts),14,16))
tpms <- cbind(tpms01A,tpms11A)
tpms_log2 <- cbind(tpms01A_log2,tpms11A_log2)

write.table(counts,"TCGA-LUAD/TCGAdata/counts.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms,"TCGA-LUAD/TCGAdata/tpms.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms_log2,"TCGA-LUAD/TCGAdata/tpms_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
```

## Day 3: 免疫得分计算与生存分析

### 几种得分计算

计算患者免疫得分与肿瘤纯度（肿瘤组分在基质+免疫+肿瘤组分中占比）
```R
library(estimate)
# Load the data; we only use the case group here so '01A'.
exp <- read.table("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 根据计算免疫得分会用到的基因与我们输入文件中现有的基因取交集进行filter
filterCommonGenes(input.f = "~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/tpms01A_log2.txt",
                  output.f = "~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/tpms01A_log2.gct",
                  id = "GeneSymbol")
# Output:
# [1] "Merged dataset includes 9883 genes (529 mismatched)."

# 全封装好了。其实我们做的主要只是整理数据
estimateScore("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/tpms01A_log2.gct",
              "~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/tpms01A_log2_estimate_score.txt",   
              platform="affymetrix")
# Output:
# [1] "1 gene set: StromalSignature  overlap= 136"
# [1] "2 gene set: ImmuneSignature  overlap= 140"


# 继续整理数据
ESTIMATE_result <- read.table("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/tpms01A_log2_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ESTIMATE_result <- ESTIMATE_result[,-1]   
colnames(ESTIMATE_result) <- ESTIMATE_result[1,]   
ESTIMATE_result <- as.data.frame(t(ESTIMATE_result[-1,]))
rownames(ESTIMATE_result) <- colnames(exp) # 发现转置完把原来的行名给变了（杠变成点了），于是把原始数据的列名赋过来

# 保存结果，就完辽
write.table(ESTIMATE_result, file = "~/Desktop/TCGAtrain/TCGA-LUAD/ESTIMATE_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F) 

```

### KM生存分析

**根据一个变量将样本划分为两组（按中位数分为高低组），看肿瘤患者的生存长短与这一变量有没有关系**。我们这个例子关心的变量是上述计算出来的免疫得分。

有两种方法可以得到生存信息：既可以从TCGA下载（没有整理好；不推荐），也可以从[XENA](https://xenabrowser.net/)的TCGA Hub下载，找到TCGA-LUAD进入phenotype - Curated survival data，提供.txt格式文件，可以直接全选并复制过来。

我们只要OS那两列。OS: overall survival，1指死亡，0指生存。
```R
survival <- read.table("~/Desktop/TCGAtrain/TCGA-LUAD/Survival_data/OS.txt", header = TRUE, sep = "\t")
# 我们只要OS那两列
survival = survival[,c(1,3,4)]
# 发现样本Barcode后缀为"01"这样子而没有最后那个"A"，因为一会儿需要与之前整理的tpms表达谱进行合并，这里也修改成"01A"形式以匹配（由于刚才导入没带参数orz，一开始读文件就没有把样本名那行当作行名，所以这里省掉rownames_to_column一步）
survival$name <- paste0(survival$sample,'A')
# 瞅一下样本名第14-16位是不是想要的形式了
table(substring(survival$name,14,16))
# Output: 
# 01A 02A 11A 
# 519   2 120
rownames(survival) = survival$name
survival = survival[,2:3]

# 合并01A生存信息与相应tpms表达谱（样本数不同，没关系，根据样本Barcode取交集。但注意这俩df一个是样本作为行另一个是作为列，所以不能直接inner_join...只能这样子生成一个common string vector再提取）
tpms01A_log2 <- read.table("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/tpms01A_log2.txt", sep = "\t",row.names = 1,check.names = F,header = T)
a <- intersect(colnames(tpms01A_log2),rownames(survival))
table(substr(a,14,16)) 
# Output:
# 01A
# 513
exp_01A <- tpms01A_log2[,a]
surv_01A <- survival[a,]
exp_01A <- exp_01A %>% t() %>% as.data.frame()
identical(rownames(exp_01A),rownames(surv_01A)) # Output: TRUE
exp_surv_01A <- cbind(surv_01A,exp_01A)
# 保存保存喵
setwd("TCGA-LUAD/TCGAdata")
write.table(surv_01A,"surv_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(exp_surv_01A,"exp_surv_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
```

我们把之前算完免疫得分那个ESTIMATE_result的数据合并过来。算得分这个也是基于tpms那个已经整理过的数据做的，很干净可以直接合并～～
```R
identical(rownames(ESTIMATE_result),rownames(surv_01A)) # Output: TRUE
ESTIMATE_result_surv_01A = cbind(surv_01A, ESTIMATE_result)
write.table(surv_01A,"~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/ESTIMATE_result_surv_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
```

然后就开始做生存分析。如上所述，先按照中位数分成免疫得分高低组：
```R
# 将生存时间的以天为单位转换为以年为单位
ESTIMATE_result_surv_01A$OS.time <- ESTIMATE_result_surv_01A$OS.time/365
# 分组：往`group`这个新建的列里填充分组标签，再转换为用于分类变量的factor
ESTIMATE_result_surv_01A$group = ifelse(ESTIMATE_result_surv_01A$ImmuneScore > median(ESTIMATE_result_surv_01A$ImmuneScore), "High", "Low")
ESTIMATE_result_surv_01A$group = factor(ESTIMATE_result_surv_01A$group, levels = c("Low", "High")) # (Recall, the **reference level** is the first level specified in the factor...)

library(survival)
fitd = survdiff(Surv(OS.time, OS) ~ group, # Surv(): creates a `survival` object
                data = ESTIMATE_result_surv_01A,
                na.action = na.exclude) 
pvalue = 1- pchisq(fitd$chisq,length(fitd$n)-1)
# 拟合生存曲线
fit = survfit(Surv(OS.time, OS) ~ group,
              data = ESTIMATE_result_surv_01A)
summary(fit)
p.label = paste0("P", ifelse(pvalue<0.001,"< 0.001", paste0(" = ",round(pvalue,3))))

# 画图
library(survminer)
ggsurvplot(fit,
           data = ESTIMATE_result_surv_01A,
           pval = p.label,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jama", # 配色，各种期刊名
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0,20), # x轴长度
           break.time.by = 5, # x轴步长
           legend.title = "ImmuneScore",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # y轴标签
           xlab = "Time (Year)", # x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()

```

![weird survival curve](https://freeimghost.net/images/2025/05/05/survival_curve.png)

完了完了P值跑出来怎么0.319比视频里面的高这么多，是哪一步数据整理错了吗完了完了完了。还是说过了两年数据库更新了...? T_T



## Day 4: 临床数据整理及绘图

还是用之前访问到的`SummarizedExperiment`对象。之前我们提取counts和TPM是从`assays`中提取了它的expression data，现在需要临床数据所以提`colData`（即phenotype data, sample covariates）。同理，S4对象用`@`提取。

```R
load("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/luad.gdc_2022.rda")
# 提取并去重
clinical <- as.data.frame(expquery2@colData) %>%   
  .[!duplicated(.$sample),]
# 提取需要的临床信息列
clinical <- clinical[,c("gender","age_at_index","ajcc_pathologic_stage","ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m")]

# 整理分类名 其实可以直接用regex：eg. "[ab]"
clinical$ajcc_pathologic_t <- gsub("a","",clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_t <- gsub("b","",clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_m <- gsub("a","",clinical$ajcc_pathologic_m)
clinical$ajcc_pathologic_m <- gsub("b","",clinical$ajcc_pathologic_m)
clinical$ajcc_pathologic_stage <- gsub("A","",clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_stage <- gsub("B","",clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_stage <- gsub("C","",clinical$ajcc_pathologic_stage)

# 整理行名
rownames(clinical) <- substring(rownames(clinical),1,16)
# 只提取01A样本（因为研究的只是肿瘤进展），与tpms表达谱合并
exp01A <- read.table("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
clinical01A <- clinical[colnames(exp01A),]
exp01A <- exp01A %>% t() %>% as.data.frame()
identical(rownames(clinical01A),rownames(exp01A)) # TRUE
clinical.expr01A <- cbind(clinical01A,exp01A)
write.table(clinical.expr01A,"clinical.expr01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# 将临床数据与ESTIMATE_result合并
ESTIMATE_result01A <- read.table("~/Desktop/TCGAtrain/TCGA-LUAD/Survival_data/ESTIMATE_result.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
identical(rownames(clinical01A),rownames(ESTIMATE_result01A)) # TRUE
clinical.ESTIMATE_result01A <- cbind(clinical01A,ESTIMATE_result01A)
# 保存为.csv，一会儿去作图
write.csv(clinical.ESTIMATE_result01A,"~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/clinical.ESTIMATE_result01A.csv")
```

之后用excel打开，选择免疫得分一列和我们关心的一列临床信息（例子中是看stage，作为解释变量）复制到空白表格中，升序排序并删掉空值行。按照StageI/II/III/IV新建四个headers，把对应stage样本的免疫得分（它就是想研究的响应变量）复制到header底下相应的列，删除原有的两列并保存。

进入[仙桃学术](https://www.xiantaozi.com/)，选择分组比较图。

![ImmuneScore stages](https://freeimghost.net/images/2025/05/05/immuneScore_stages.png)

// 太好了终于跟视频里演示的一样了。Q_Q



## Day 5: 差异分析

每次都把数据整理好，之后用的时候索引名及排列顺序都是一样的，就可以很方便～～

还是用免疫得分来分高低组

```R
library(DESeq2)
library(tidyverse)
setwd("TCGA-LUAD/TCGAdata")
counts_01A <- read.table("counts01A.txt", sep = "\t", row.names = 1, check.names = F, header = T)
estimate <- read.table("ESTIMATE_result.txt", sep = "\t", row.names = 1, check.names = F, header = T)

# 分组
x = "ImmuneScore" # 这里换成'StromalScore'等等estimate文件中的其他列名不就可以根据其他相应得分进行样本分组了嘛
med = as.numeric(median(estimate[,x]))
estimate = as.data.frame(t(estimate))
identical(colnames(counts_01A),colnames(estimate))
conditions=data.frame(sample=colnames(counts_01A), group=factor(ifelse(estimate[x,]>med,"high","low"),levels = c("low","high"))) %>% column_to_rownames("sample")

# 用DESeq2跑差异分析
if (!is.matrix(counts_01A)) counts_01A <- as.matrix(counts_01A)
dds <- DESeqDataSetFromMatrix(
  countData = counts_01A,  # counts的作用就只在这里
  colData = conditions,
  design = ~ group)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
save(res,file = "DEG_ImmuneScore.rda")
```

#### 理解、定义并筛选差异基因

看logFC值正负时一定要整明白**讲上/下调是相对于什么条件**。

回忆上面创建`conditions`这个factor作为分类变量时是设置了`levels = c("low","high")`，因为还是按顺序决定所以low组（排在第一个的level）自动作为reference level，相当于对照组。

- **“上调”就是指表达量高于对照组，即在high组中表达显著升高**。

- 下调同理，eg. `logFC = -1`表示high组中该基因的表达量是low组的2^(-1)=0.5倍。

注意注意，我们下面筛选差异基因的步骤中会只保留logFC绝对值大于1的（就是跟对照组比表达量为两倍/一半）。这样子做其实有可能漏掉一些基因，它们或许表达差异不到两倍但影响很大。后面在Day8我们做GSEA的时候就不这样子筛选了，用全部DEG去跑。

```R
# 定义并筛选差异基因
exp <- read.table("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/tpms01A_log2.txt", sep = "\t", row.names = 1, check.names = F, header = T)
DEG <- as.data.frame(res)
logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT")) # 新建一列，填入差异表达信息
table(DEG$change) 
# DOWN   NOT    UP 
#  553 18230   708 

a <- filter(DEG,change == 'UP')
b <- filter(DEG,change == 'DOWN')
c <- rbind(a,b)
exp_diff <- exp[rownames(c),]
```

#### 热图绘制及表达谱样本整理

注意此时**表达矩阵的样本/列之间**没有按照组别（免疫得分高/低组）进行整理，所以现在看样本间是乱序的。
```R
# 用筛选得到DEG的表达谱绘制热图
library(pheatmap)
pheatmap(exp_diff,
         annotation_col = conditions, # 图上方那个表示样本的彩色横框
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows = F,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)

```
![unsortedDEGheatmap](https://freeimghost.net/images/2025/05/05/unsortedDEGheatmap.png)


对表达谱的列根据样本分组进行整理，绘出的热图就能看到**在不同条件下差异基因的表达上下调情况**。再加上按行聚类。

```R
exp_diff_high = exp_diff[,conditions$group=='high']
exp_diff_low = exp_diff[,conditions$group=='low']
exp_diff_sorted <- cbind(exp_diff_high,exp_diff_low)

pheatmap(exp_diff_sorted,
         annotation_col = conditions, # 注意看这个！现在是整好了的
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
```
![sortedDEGheatmap](https://freeimghost.net/images/2025/05/05/sortedDEGheatmap.png)
（太丑了，dbq，不应该加聚类的55555 有点可怕了）


### 两次差异分析结果取交集

需要先以其它条件再做一次差异分析。这次不用免疫得分了，用基质得分划分高低组（作为不同样本条件）来寻找差异基因。过程完全照搬：

```R
library(DESeq2)
library(tidyverse)
library(pheatmap)
setwd("TCGA-LUAD/TCGAdata")
counts_01A <- read.table("counts01A.txt", sep = "\t", row.names = 1, check.names = F, header = T)
estimate <- read.table("ESTIMATE_result.txt", sep = "\t", row.names = 1, check.names = F, header = T)
# 给样本分组
x = "StromalScore" 
med = as.numeric(median(estimate[,x]))
estimate = as.data.frame(t(estimate))
identical(colnames(counts_01A),colnames(estimate))
conditions=data.frame(sample=colnames(counts_01A), group=factor(ifelse(estimate[x,]>med,"high","low"),levels = c("low","high"))) %>% column_to_rownames("sample")
# DESeq跑差异分析
if (!is.matrix(counts_01A)) counts_01A <- as.matrix(counts_01A)
dds <- DESeqDataSetFromMatrix(
  countData = counts_01A,
  colData = conditions,
  design = ~ group)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
save(res,file = "DEG_StromalScore.rda")
# 筛选差异基因
exp <- read.table("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/tpms01A_log2.txt", sep = "\t", row.names = 1, check.names = F, header = T)
DEG <- as.data.frame(res)
logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT")) # 
table(DEG$change) 
# Output:
# DOWN   NOT    UP 
# 496 18394   601 
a <- filter(DEG,change == 'UP')
b <- filter(DEG,change == 'DOWN')
c <- rbind(a,b)
exp_diff <- exp[rownames(c),]
exp_diff_high = exp_diff[,conditions$group=='high']
exp_diff_low = exp_diff[,conditions$group=='low']
exp_diff_sorted <- cbind(exp_diff_high,exp_diff_low)
# 绘图
pheatmap(exp_diff_sorted,
         annotation_col = conditions,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
```
![StromalDEGheatmap](https://freeimghost.net/images/2025/05/05/StromalDEGheatmap.png)

现在来取交集。先把两次差异分析跑出来的这些上下调的基因提取出来。
```R
# 上述代码不是已经跑过这个了嘛，环境变量也都还在
# a <- filter(DEG,change == 'UP')
# b <- filter(DEG,change == 'DOWN')
# 所以直接写入
write.csv(a, file = "~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/Stromal_up.csv")
write.csv(b, file = "~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/Stromal_down.csv")

# 按ImmuneScore分组的同理

```
在excel中打开这些csv，选中两组差异分析中得到的上/下调基因名，新建两个表格分别放上调与下调，按列复制过去（所以一共是两列，一个Stromal一个Immune），保存。

再打开仙桃学术，选择韦恩图，除绘图外**它可以输出"交集情况.xlsx"，导出**。

![Venn for DEG_down](https://freeimghost.net/images/2025/05/06/Venn-for-DEG.png)

excel中把仙桃学术给上下调基因分别取到的交集（两个"交集情况.xlsw"文件的第三列）贴在同一列，将这个单独一列保存为txt，导入r。

### 用交集差异基因做GO和KEGG富集分析

> **"Are genes in the predefined set more frequently observed in my results than expected by random chance?"** -> quantified by comparing the observed overlap between the sets to a null distribution.

回忆，我们用*按照肿瘤患者的免疫得分和基质得分来分组做的两次差异分析*得到的差异基因又取了交集，相当于更进一步筛选相关基因。然后我们就做富集分析，GO和KEGG两个都用一下（都需要联网去数据库上fetch），后者富集的是通路。哇生竞记忆又回来了-_-，当时那堆生信的怪名词怪缩写大概确实没白背，虽然现在还是会觉得how upsetting。

```R
union <- read.csv("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/union.txt", sep="")

# 回到我们的DESeqResult对象（用哪次差异分析得到的都可以），在里面提取出交集差异基因的相关数据。
DEG = as.data.frame(res)
DEG = DEG[union$SYMBOL,] 
DEG =  rownames_to_column(DEG, 'SYMBOL')
library(org.Hs.eg.db) # 又是它，'genome wide annotation for Human, primarily based on **mapping** using Entrez Gene identifiers.' 注意注意！！GO和KEGG只能识别Entrez ID不能识别SYMBOL这种形式，所以我们才要在这里map一下。

library(clusterProfiler) # 这个就是做富集分析的包捏
# `bitr() function from `clusterProfiler` package: "Biological Id TRanslator"
genelist = bitr(DEG$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") 
DEG = inner_join(DEG,genelist,by='SYMBOL')

# GO分析
ego = enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all", # GO中的那三类(MF/CC/BP)都要
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)
ego_res <- ego@result
save(ego,ego_res,file = "~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/GO_DEG_union.Rda")

kk = enrichKEGG(gene = DEG$ENTREZID,
                organism = 'hsa',
                pvalueCutoff =0.1, 
                qvalueCutoff =0.1)
kk_res <- kk@result
save(kk,kk_res,file = "~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/KEGG_DEG_union.Rda")
```

#### 富集分析后绘制网络图

```R
library(ggnewscale)
logFC = DEG$log2FoldChange # 把表达差异这一列从差异分析结果中提取出来作为一个向量
names(logFC)= DEG$ENTREZID # 用对应基因的Entrez ID作为向量每一元素（即该基因表达差异）的索引，与富集分析结果进行匹配
logFC = sort(logFC,decreasing = T)
head(logFC)
#     5545     1178     8115    93978   441168     2208 
# 3.107605 2.616990 2.177435 2.150608 2.138932 2.117131 

cnetplot(ego, foldChange = logFC, circular = TRUE, colorEdge = TRUE, layout = "circle")
cnetplot(kk, foldChange = logFC, circular = TRUE, colorEdge = TRUE, layout = "circle")
# 画出来标的基因ID是SYMBOL的形式不是Entrez ID
```
为什么跑出来的图跟视频演示又这么不一样Q_Q

![enrichKEGG](https://freeimghost.net/images/2025/05/06/enrichKEGG.png)

![enrichGO](https://freeimghost.net/images/2025/05/06/enrichGO.png)


## Day 6: 蛋白质互作及COX回归分析

### PPI

#### STRING与Cytoscape

突然就想起来南昌集训的时候生信老师还展示过Cytoscape，当时觉得真是好高级啊。果然还是很怀念Q_Q，非常非常非常怀念。

好了请止住你的碎碎念。

[STRING](https://cn.string-db.org/cgi/input?sessionId=bXmYsv7CnUrH&input_page_active_form=multiple_identifiers)中'Search' -> 'Multiple Proteins by Names / Identifiers'，粘贴过去我们之前取交集取到的一列基因名 -> 'The following proteins in Homo sapiens appear to match your input. Please review the list, then click 'Continue' to proceed.' -> 返回一张蛋白作节点的网络图。可以修改置信度来筛选蛋白，让显示的网络中的节点不要那么多，然后**导出为.tsv格式**，可以在Excel或Cytoscape中打开并进一步整理（现在的图太丑了）。

**Cytoscape**中'Import' -> 'Network from File'导入.tsv，会生成网络图。

- Cmd方框选择，Cmd+Shift套索工具，选中节点后可以鼠标拖动。

- **'Tools' menu bar -> 'Analyze Network' 返回网络分析表格** -> 'Export Table'导出。

- 'Style'中修改网络图样式，'Fill Color'中可以根据'Column'中各种参数对节点颜色进行不同的mapping，这里选择根据degree（该节点度数，即跟它有互作的蛋白数）填充颜色。可以发现CD4度数很大，这种基因称*hub gene*。应该导出为pdf（下图的png就非常不好，但是我不再去整一遍了）T_T

![PPI network from Cytoscape](https://freeimghost.net/images/2025/05/07/PPI_network_from_Cytoscape.png)

在Excel中提取出从Cytoscape中导出的网络分析表格的Degree列和基因名一列，在新表格中根据Degree按降序排序。按照我们要复现的文章，只保留top 30的基因，进入仙桃学术绘画一维柱状图。发现x轴标注都挤到一起去了，可以调整参数中'坐标轴'的'x轴标注旋转'及在'图片'选项调整图像宽度。

![PPI top30](https://freeimghost.net/images/2025/05/07/PPI_top30.png)

### 单因素COX回归分析

**寻找基因表达量与肿瘤患者生存的关系**，尤其看生存时间，因为Cox回归（全称为Cox比例风险模型）主要探讨终点事件发生速度有关的因素，这是它与使用logistic regression这种二分类场景很不同的一点。

复现中仍然使用我们之前筛选的差异基因（两次差异分析后取了交集的），从之前整理的表达谱-生存信息数据中提取出对应基因的数据，逐个遍历计算**HR值（大于1代表该基因对患者是一个危险因素，小于1代表是保护因素）** 及相应检验统计量、显著性、HR值的95%置信区间等等，最后绘制HR森林图。

```R
library(survival)
library(forestplot)
DEG_union <- read.csv("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/union.txt", sep="")
exp_surv_01A <- read.table("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/exp_surv_01A.txt", sep = "\t", row.names = 1, check.names = F, header = T)
# 提取出患者生存信息（前两列），再用取到交集的基因名作为索引提取对应基因表达谱的列（注意读进来这俩表是转置的）
surv.expr = cbind(exp_surv_01A[,1:2],exp_surv_01A[,DEG_union$SYMBOL])

# 把每一个基因（每一列）都遍历一遍，跑单变量COX模型
Coxoutput <- NULL 
for(i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
  cox <- coxph(Surv(OS.time,OS) ~ surv.expr[,i], data = surv.expr)
  coxSummary = summary(cox)
  
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]), 
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}

# 按p值筛选基因
Coxoutput = arrange(Coxoutput,pvalue)
topgene = Coxoutput[Coxoutput$pvalue < 0.005,] # 按照复现文章设置的阈值，现在从566个filter剩下51个
write.csv(topgene, file = "~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/cox_gene_sig.csv")
```

### 森林图绘制及读图
```R
# 制作输入表格
tabletext <- cbind(c("characteristics",topgene$gene),
                   c("HR",format(round(as.numeric(topgene$HR),3),nsmall = 3)),
                   c("lower 95%CI",format(round(as.numeric(topgene$lower),3),nsmall = 3)),
                   c("upper 95%CI",format(round(as.numeric(topgene$upper),3),nsmall = 3)),
                   c("pvalue",format(round(as.numeric(topgene$p),3),nsmall = 3)))

forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(topgene$HR)),
           lower=c(NA,as.numeric(topgene$lower)), 
           upper=c(NA,as.numeric(topgene$upper)),
           graph.pos=5, # 图在表中的列位置
           graphwidth = unit(.25,"npc"), # 图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI", # box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"), # box颜色
           
           boxsize=0.4, # box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=T, # 显示区间
           zero=1, # zero线横坐标
           lwd.zero=1.5, # zero线宽
           xticks = c(0.5,1,5), # 横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1.2), # 各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"), # 在第一行上面画黑色实线
                           "2" = gpar(lwd=1.5, col="black"), # 在第一行标题行下画黑色实线
                           "53" = gpar(lwd=2, col="black")), # 在最后一行（第53行，因为一共51个基因加标题行）上画黑色实线
           lineheight = unit(.75,"cm"), # 固定行高
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
```

其实非常好读的说，就完全是可视化了一下那个输入表格～～。那条竖线（HR = 1）就是我们之前说的危险/保护因素分界线。还可以再去Adobe Illustrator把图上过了竖线的基因/box改成红色的，危险因素嘛。

![COX forestplot](https://freeimghost.net/images/2025/05/07/COX_forestplot945969790948b851.png)

#### 解决图像显示/导出不完整问题！！

直接运行绘制森林图的那行代码，在RStudio小窗口中怎么调也显示不全，导出之后还是不全。找到了一个解决办法：不要在窗口中打开，然后这样子跑
```R
png("file_name.png",height = 1300,width = 800) # 设置好保存尺寸
forestplot(...)
dev.off()
```
这样子就可以直接完整地保存到目标路径！！！！终于整出来完整的图了，很开心，必须贴出来，就在上面读图那里。


### 取PPI与Cox分别筛选后的基因交集

之前STRING做蛋白质互作的时候返回一张表，不是按照节点度数前三十筛选了基因嘛。然后上述Cox跑完了也根据p值筛选到了那五十来个过了阈值的。我们取这俩set的交集，还是一样把基因名按列粘贴到同一个新的表格中，又去仙桃学术里画韦恩图，然后导出生成的"交集情况.xlsx"。取到的交集有五个基因，可怕的是发现要复现的文章下文做的单基因分析用到的BTK基因在我的交集里面没取到，Cox筛选到它了但PPI前三十的没有它，啊啊Q_Q（文章里取到的交集有两个基因，BTK与CCR2，后者这次倒是也取到了...）

![COX_PPI_Venn](https://freeimghost.net/images/2025/05/07/COX_PPI_Venn.png)

## Day 7-8: 用筛选到的单基因BTK做各种分析

如上，这些单基因分析都是用BTK基因做的q_q，虽然没在复现取到的交集里面但就这样跟着做吧5555。

### 正常/肿瘤样本中表达量比较-柱状图

分别读入之前整理好的01A与11A样本的tpms表达谱数据，通过行索引`["BTK",]`提取出该基因在全部样本中表达量的一整行，再转置成列。把这两列（一列是01A的全部样本一个是11A的）拼到一起，不要样本Barcode信息只留表达量的数值，两列的headers分别改成'Normal'和'Tumor'，这样子很快速就整理完了。

再将这一表格导入仙桃学术画分组比较图中的柱状图，两个柱子就对应两列（'Normal'和'Tumor'）。参数调一调，导出。

![Lecture screenshot](https://freeimghost.net/images/2025/05/17/Screenshot-2025-05-17-at-16.18.00.png)

### 样本配对原理及配对图

**能配对是因为同一患者的肿瘤标本周围会带正常组织，取标本的时候有时会把肿瘤患者的正常组织也取上并送去测序，所以11A中是有来自同一样本的故可以配对上的**（一会儿提取样本的时候就能发现11A的58行中有57行都是在交集里面的）。匹配样本Barcode（除了01A/11A字符）就可以。其他的整理步骤跟上述柱状图很类似

```R
tpms01A_log2 = read.table("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tpms11A_log2 = read.table("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/tpms11A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tpms01A_log2 <- tpms01A_log2 %>% t() %>% as.data.frame()
tpms11A_log2 <- tpms11A_log2 %>% t() %>% as.data.frame()
rownames(tpms01A_log2) <- substring(rownames(tpms01A_log2),1,12)
rownames(tpms11A_log2) <- substring(rownames(tpms11A_log2),1,12)
a <- intersect(rownames(tpms01A_log2),rownames(tpms11A_log2))
tpms01A_log2 <- tpms01A_log2[a,] 
tpms11A_log2 <- tpms11A_log2[a,] # 这两行提取用的行索引是一样的，所以接下来拼完了顺序也是配对的。
paired <- cbind(tpms11A_log2$"BTK",tpms01A_log2$"BTK") # 11A放在前面
paired <- as.data.frame(paired)
write.csv(paired,file = "paired_BTK.csv")
```

接着拿到表格之后还是像柱状图那样在excel里改两个列的headers，然后导入到仙桃学术里面绘制配对图，y轴标注忘改了，orz。

![paired_BTK](https://freeimghost.net/images/2025/05/07/paired_BTK.png)

### 用单基因表达量分高低组做生存分析

把之前用免疫得分分组的代码复制过来，改一改：
```R
surv <- read.table("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/exp_surv_01A.txt", sep = "\t", row.names = 1, check.names = F, header = T)
surv$OS.time <- surv$OS.time/365
surv$group = ifelse(surv$BTK > median(surv$BTK), "High", "Low") # 主要改动就分组这里
surv$group = factor(surv$group, levels = c("Low", "High")) 
table(surv$group) 
# Low High
# 257 256
library(survival)
fitd = survdiff(Surv(OS.time, OS) ~ group,
                data = surv,
                na.action = na.exclude) 
pvalue = 1- pchisq(fitd$chisq,length(fitd$n)-1)
fit = survfit(Surv(OS.time, OS) ~ group,
              data = surv)
p.label = paste0("P", ifelse(pvalue<0.001,"< 0.001", paste0(" = ",round(pvalue,3))))
library(survminer)
ggsurvplot(fit,
           data = surv,
           pval = p.label,
           conf.int = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata",
           palette = "jama", 
           legend.labs = c("Low", "High"), 
           size = 1,
           xlim = c(0,20),
           break.time.by = 5, 
           legend.title = "BTK Expression",
           surv.median.line = "hv", 
           ylab = "Survival probability (%)", 
           xlab = "Time (Year)", 
           ncensor.plot = TRUE, 
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
# dev.off()
```

![BTK survival curve](https://freeimghost.net/images/2025/05/07/BTK_survival_curve.png)

这个画出来p值终于小了！！！！

### 比较不同临床分期的样本中BTK表达量

找到我们之前整理好的临床信息-表达量数据，提取出全部样本的临床信息（在第1至6列）及BTK（这个整理好的数据表中基因也作列！）在各样本中的表达量，保存为新表格。
```R
clinical.expr01A <- read.table("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/clinical.expr01A.txt", sep = "\t", row.names = 1, check.names = F, header = T)
clinical_BTK = cbind(clinical.expr01A[,1:6],clinical.expr01A$BTK)
write.csv(clinical_BTK, file = "~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/clinical_BTK.csv")
```

然后又是去excel里面把想比较的列（想研究的变量）提出来放到新表格里，按照不同分期分成几列，反正整理和之后去仙桃学术绘制分组比较图的流程与上文Day4那里完全相同～～

### BTK表达量高低组的差异分析及富集分析

做差异分析跟上文用免疫得分高低划分样本组的完全一样，只是分组的条件判断变成了跟BTK表达量的中位数比（同刚才跑的生存分析）。

我的DESeq2包中`DESeqDataSetFromMatrix()`这两天用不了了，明明前些天还在学coursera的时候还好好的呢Q_Q，跑那会儿一样的代码也不行。昨天早上试了什么办法都还是报错，后来在Bioconductor的support上面看到有人在几小时前也提了同样的问题，大概是新的bug，先等等吧5555，反正流程是完全一样的就偷一下懒直接load进来老师跑好的：`load("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/DEG_BTK.Rda")`，于是环境变量里得到`DESeqResults`对象'res'一只。

#### 用到本地基因集gmt注释文件的GSEA

差异基因拿到了接着跑基因组富集分析，多数代码还是跟Day5很像。

我们这次GSEA与上次用GO和KEGG的另一个不同是！！它需要下载到本地一个gmt文件（从MSigDB (Molecular Signatures Database)数据库里下载的注释，[整个文件夹将所有基因集划分成八类](https://zhuanlan.zhihu.com/p/504101161)），这样子就不像上次的两个跑的时候必须联网到相应的数据库里fetch到基因集，而且也不用因为GO和KEGG只识别Entrez ID不识别SYMBOL形式进行一步用到了org.Hs.eg.db包的基因名映射（这一步还有可能会丢掉一些不能转换的基因），因为MSigDB文件里面已经有SYMBOL版的注释了。不信的话可以打开文件夹看看www
```R
load("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/DEG_BTK.Rda")
DEG = as.data.frame(res) %>% arrange(padj) %>% # 19934个差异基因，按照p值从小往大排
  rownames_to_column('Gene') # 这个'Gene'列底下放的还是SYMBOL形式的基因名喵

library(clusterProfiler)

geneList = DEG$log2FoldChange
names(geneList) = as.character(DEG$Gene)
geneList = sort(geneList, decreasing = TRUE)
head(geneList)
# PRB4    RXFP2      CLC  CEACAM8     AZU1   MCEMP1 
# 3.120753 2.672972 2.493237 2.244911 2.058411 2.058164 

# 虽然下载了整个MSigDB的文件夹下来，其实只用其中的这一个gmt文件，但是也学一下这样子指定路径的gmt文件如何读取！！
msigdb_GMTs <- "msigdb_v7.0_GMTs" # 这个文件夹就放在project路径里面
msigdb <- "h.all.v7.0.symbols.gmt" # 这样子指定路径以后想用其他注释文件的时候方便更改呢！
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))
# 仅两列，每一行'term'字段说这是什么通路，'gene'字段放相应SYMBOL基因名
head(table(kegmt$term)) # 每个通路出现了多少行，即这个一样的通路中有多少基因
#   HALLMARK_TNFA_SIGNALING_VIA_NFKB 
#                                200 
#                   HALLMARK_HYPOXIA 
#                                200 
#   HALLMARK_CHOLESTEROL_HOMEOSTASIS 
#                                 74 
#           HALLMARK_MITOTIC_SPINDLE 
#                                199 
# HALLMARK_WNT_BETA_CATENIN_SIGNALING 
#                                 42 
#        HALLMARK_TGF_BETA_SIGNALING 
#                                 54 

set.seed(1) 
gsea <- GSEA(geneList,TERM2GENE = kegmt) # 跑GSEA
gsea_result_df <- as.data.frame(gsea)
nrow(gsea_result_df) # 29，即从我们的19934个差异基因找到了29条统计学显著的通路。
save(gsea,gsea_result_df,file = "GSEA_BTK_h.all.rda") # 保存为Rdata很好用的一点就是，可以把几个变量一块保存进去


# 上面是用HALLMARK collection去找的基因集。接着再用C7 collection，就是免疫相关基因集那个。如上所述，gmt的指定路径更改一下就可以力其他代码不用变呢
msigdb2 <- "c7.all.v7.0.symbols.gmt" # 原来那个注释文件叫h.all.v7.0.symbols.gmt来着
kegmt2 <- read.gmt(file.path(msigdb_GMTs,msigdb2))
set.seed(1) 
gsea2 <- GSEA(geneList,TERM2GENE = kegmt2) 
gsea_result_df2 <- as.data.frame(gsea2)
nrow(gsea_result_df2) # 1791...??? 比Hallmark中富集到的29个多了好多好多
save(gsea2,gsea_result_df2,file = "GSEA_BTK_c7.all.rda") 
```

#### 富集分析绘图并读图

认为|NES|>1，p-val<0.05，FDR q-val<0.25是显著富集的，其中NES为归一化后的ES值。

```R
# 单结果绘制
library(enrichplot)
gseaplot2(gsea,8,color="red",pvalue_table = T) # 画出来第八行的通路，可以从gsea_result_df中看到该行ES值为正所以曲线是凸的
gseaplot2(gsea,9,color="red",pvalue_table = T) # 第九行ES值就为负，所以曲线concavity也相反

# 多个结果绘制
gseaplot2(gsea, geneSetID = c(8,10,12), subplots = 1:3)
gseaplot2(gsea, geneSetID = 1:3, subplots = 1:3)

dev.off()
```

![GSEA enrichplot](https://freeimghost.net/images/2025/05/07/GSEA_enrichplot.png)

![negative_ES enrichplot](https://freeimghost.net/images/2025/05/07/negativeNES_enrichplot.png)

**如何读图**：

- **富集曲线反映该基因集/通路中基因在排序列表中的分布**，峰越高 ('Running Enrichment Score'越大），说明富集程度越强。从左至右每到一个基因（**注意我们的geneList是按照logFC排的**），计算出一个ES值，连成线。底下那些竖线就是每次geneList中一个基因在该条通路中*ES为正就打一条*，所以竖线分布越密集也就对应此通路在整个geneList中富集的位置。

  - ES值为正：曲线整体趋势从左侧上升。表示基因集中的基因 **富集在logFC较大的区域（排在geneList前边），即该通路中的基因在处理组条件下显著上调（正调控）！！** 一定回忆差异分析的时候说logFC正负和对照组/处理组怎么看

  - ES值为负：从左侧下降（呈V型这样子），表示该通路/基因组中的基因在处理组条件下显著下调（负调控），或者反过来说就是在对照组（这里是BTK表达量低组）中富集。
  
  - 底部热图/颜色条是geneList**排序指标**，也就是按照logFC嘛，红色为高logFC（对照组条件下的上调基因），蓝色为低/负logFC（下调基因）。


**发现BTK高组中富集到了一些免疫通路，低组富集到了一些代谢通路，这就预示着肿瘤微环境在进展过程中发生了一些变化。**

## Day 9: CIBERSORT免疫浸润分析

计算每一个样本中22种免疫细胞各占的百分比，需要一个R文件存放CIBERSORT的算法和一个'LM22'计算各免疫细胞的参考文件：
```
## Gene symbol  B cells naive B cells memory Plasma cells T cells CD8 T cells CD4 naive
## ABCB4        555.71          10.74        7.226       4.311             4.606
## ABCB9         15.60          22.09      653.392      24.224            35.672
## ACAP1        215.31         321.62       38.617    1055.613          1790.097
## ACHE          15.12          16.65       22.124      13.428            27.188
## ACP5         605.90        1935.20     1120.105     306.313           744.657
```

```R
library(e1071)
library(parallel)
library(preprocessCore)
library(tidyverse)
source("CIBERSORT.R")   # 这样子导入R文件
sig_matrix <- "LM22.txt"   
mixture_file = '~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/tpms01A_log2.txt'
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm = 100, QN=TRUE) # `perm = 100`: 计算100次，次数越多结果越收敛
res_cibersort <- res_cibersort[,1:22]  # 只要前22列关于22种免疫细胞的信息，后面三列的'P-value', 'Correlation', 'RMSE'就不要了
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   #  去除丰度在所有样本（即这一整列的和）全为0的细胞，filter完还是22列
ciber.res <- as.data.frame(ciber.res)
write.table(ciber.res,"~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/ciber.res.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
```

### 彩虹图绘制及读图

![彩虹图](https://freeimghost.net/images/2025/05/07/CIBERSORT_rainbowplot.png)

如何读图看下面的代码注释，其实**只是把相对比例可视化出来了**（但注意它是**样本作列**细胞作行，因为把之前得到的ciber.res转置了）。

```R
mycolor <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) # 按照ciber.res的列数（即免疫细胞种类数）创建彩虹色板，带70%透明度
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
# 一根一根的柱子，一共513列 **每一列代表一个样本**（df需要先转置）：
barplot(as.matrix(t(ciber.res)),
        border = NA, # 柱子无边框，所以像这样子每根彩虹的柱子stack起来
        names.arg = rep("",nrow(ciber.res)), # 无横坐标样本名
        yaxt = "n", # 先不绘制y轴，下一行代码axis()那个补充y轴标注
        ylab = "Relative percentage", # 我们ciber.res算出来的是22种免疫细胞的相对比例嘛！！有没有感觉图很好看懂了
        col = mycolor) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%")) 
legend(par("usr")[2]-20,   
       par("usr")[4], 
       legend = colnames(ciber.res), # 添加图例，即各个颜色对应哪种免疫细胞
       xpd = T,
       fill = mycol,
       cex = 0.6, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
```

### 相关性热图

其实可视化想展示出来的内容跟上面那个彩虹图是一样的。

虽然我们这里画的相关性热图是同一个矩阵跟自己比（**看各列/各变量之间相关性**，现在这个ciber.res行列很像是常规数据分析中那种design matrix就是样本作行特征作列），但是需要回顾的是**算相关性矩阵的两个输入矩阵每一行行所代表的样本必须是一模一样的**。

```R
library(corrplot)
library(ggcorrplot)
# 注意看！！sapply中第二个参数这种情况该怎么传：
cor <-sapply(ciber.res,function(x,y) cor(x,y,method="spearman"),ciber.res)
rownames(cor) = colnames(ciber.res) # cor这个相关性矩阵原本没行索引，反正原始矩阵是自己跟自己的列比，所以行名就是原矩阵列名（细胞类型）

png("./plots/correlation_matrix_visualization.png",height = 700,width = 700)
ggcorrplot(cor, 
           hc.order = TRUE, # 使用hc.order进行排序
           type = "lower", # 热图位置在右下角还是左上角，现在输出的就在右下角
           outline.color = "white", # 轮廓颜色
           lab = TRUE, # 在热图上面标上相关系数那堆数值
           ggtheme = ggplot2::theme_gray, # 默认值为theme_minimal
           colors = c("#01468b", "white", "#ee0000"))
dev.off()
```

![correlation matrix visualization](https://freeimghost.net/images/2025/05/08/correlation_matrix_visualization.png)



### BTK表达量高低组的免疫细胞分组比较图

前面数据整理环节还是这些步骤～～

但是但是但是！在`ggpubr::ggboxplot()`画图前的数据整理还挺不一样的：注意`gather(a,key=CIBERSORT,value = Fraction,-c(group,sample))`这行的`gather()`函数！！！

-  Reshapes the data from wide format to **long format (tidy data) for plotting** (513 rows -> 513 * 22 = 11286 rows; 24 columns -> 4 columns). Takes all columns except group and sample (`-c(group, sample)`) into two new columns: `CIBERSORT` to store the immune cell type names and 
`Fraction` to store the corresponding cell fraction value for each immune cell type. 

  - Output `b`: a long-format dataframe with 4 columns:
```  
> head(b)
             sample group     CIBERSORT    Fraction
1 TCGA-44-8120-01A   low B cells naive 0.030837181
2 TCGA-99-8025-01A   low B cells naive 0.033631422
3 TCGA-50-6594-01A   low B cells naive 0.000000000
4 TCGA-78-7161-01A   low B cells naive 0.001200799
5 TCGA-55-7903-01A   low B cells naive 0.065245199
6 TCGA-38-4632-01A  high B cells naive 0.024616157
```

想一想，毕竟要作为绘图函数的输入数据，得把**想研究的变量**用**输入数据的列**表示嘛。

```R
exp <- read.table("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
med=median(as.numeric(exp["BTK",]))
exp <- exp %>% t() %>% as.data.frame()
exp <- exp %>% mutate(group=factor(ifelse(exp$BTK>med,"high","low"),levels = c("low","high"))) # Recap: "mutate() **creates new columns** that are functions of existing variables. It can also modify (if the name is the same as an existing column) and delete columns (by setting their value to NULL)."
class(exp$group) # "factor"
a <- ciber.res
identical(rownames(a),rownames(exp)) # TRUE
a$group <- exp$group
a <- a %>% rownames_to_column("sample")
library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(a,key=CIBERSORT,value = Fraction,-c(group,sample))
png("./plots/ggboxplot.png",height = 700,width = 700)
ggboxplot(b, x = "CIBERSORT", y = "Fraction", # 变量一定跟b中的相应**列**对应起来
          fill = "group", palette = "jco")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()
```
![ggboxplot](https://freeimghost.net/images/2025/05/08/ggboxplot.png)

可以看一下哪些细胞类型在高低组的**差异**是统计学显著的


### BTK单基因表达量与免疫细胞占比 相关性散点图

注：关于读入表格数据行/列索引带空格的问题：

- r默认不支持索引带空格，但读取数据时有时候会自动替换空格为`.`或`_`，然后就可以这样子引用：`B.cells.naive`（可以先eg. `colnames(df)`瞅一下）

- 列名可能保留空格，这时必须用**反引号**引用！！见下。

```R
exp <- read.table("~/Desktop/TCGAtrain/TCGA-LUAD/TCGAdata/tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) # 上个代码块已经读取到环境里了，但是整理过程中把原始数据覆盖掉了啊啊。我觉着以后还是不要动不动就把原始数据覆盖掉了...? 有的读进来挺花时间的
exp = exp['BTK',] # 只取出BTK的表达量，顺便再啰嗦一句生信里面这些expression matrix之所以会让样本作列不作行（所以我们现在取基因取的是行）是因为测序的基因数远超过样本量对吧，而r处理行的能力远超过列！！
ciber <- ciber.res %>% t() %>% as.data.frame() #转置一下，因为记得不ciber.res是更常规的design matrix那种样本作行的形式。 赋给一个新变量...这次不直接覆盖掉了
identical(colnames(ciber),colnames(exp)) # TRUE
exp_ciber <- rbind(exp,ciber)
exp_ciber <- exp_ciber %>% t() %>% as.data.frame() # 再转置，又变回design matrix形式；还是这句话，绘图中关心的变量得是输入数据的列嘛！
library(ggstatsplot)
library(ggside)
png("./plots/Naive_B_cell_correlation_scatterplot.png",height = 700,width = 700)
ggscatterstats(data = exp_ciber, # 要分析的数据
               y = BTK, # 设置y轴，为BTK表达量这一列
               x = `B cells naive`, # 设置x轴，比如这里我们先看BTK表达量与初始B细胞的关系。注意引用方式！！本来设置的变量x直接是B.cells.naive，报错了说没找到。
               type = "nonparametric", 
               margins = "both", # 是否显示边缘，默认为true                                
               xfill = "#01468b", # x轴边缘图形的颜色
               yfill = "#ee0000", # y轴边缘图形的颜色
               marginal.type = "densigram") # 在图片坐标轴边缘添加图形类型
dev.off()
```

![Naive B cell correlation scatterplot](https://freeimghost.net/images/2025/05/08/Naive_B_cell_correlation_scatterplot.png)

哇！！！刚好还是个负相关的呢！！

#### 一定区分上述的差异分析与相关性分析

刚才这个相关性分析图画出来也会返回一个p值，所以我们也还可以研究哪些细胞类型与BTK表达量的**相关性**是显著的，一定跟上个"BTK表达量高低组的免疫细胞分组比较图"返回的**差异**显著的细胞类型区分开啊（类似于independence test和homogeneity test很不同，反正区分）！

这两个分析 **（一个差异分析一个相关性分析）** 能分别筛选出不同的两组免疫细胞，文章中也给这些细胞取了交集画了韦恩图，这里就不复现那张辽（溜走）。

---

## 总结：文章讲了什么故事...?

**文章标题**：BTK Has Potential to Be a Prognostic Factor for Lung Adenocarcinoma and an Indicator for Tumor Microenvironment Remodeling: A Study Based on TCGA Data Mining 

（虽然旧了而且比较水，但是思路对于小白很值得学习）

从Day3到Day9复现了这么多图，所以文章到底想讲一个什么样的故事？

- 先是用免疫/基质得分做生存分析、比较临床分期，发现它们跟肿瘤预后有显著相关性，所以**确定免疫/基质得分作为入手点来筛选基因**。

- 于是就用两个得分给样本划分高低两组做差异分析找到了不同的两组差异基因，之后又去了交集。用交集中筛选出来的差异基因做了GO和KEGG分析。

- 又用这些差异基因做了蛋白质互作分析，取了top30基因。还是用交集里的这些基因，做COX单因素回归分析，取出和预后相关p值小于阈值的51个基因，跟PPI top 30的取了交集，这样双重筛选下来文章筛选到了两个基因，其中一个是BTK基因，**这就是文章怎么筛到了BTK这个基因**。

- **筛选到了单基因之后就是BTK的单基因分析了**：看正常组织/肿瘤组织中表达差异、生存分析和比较临床分期等**看它的表达量跟预后的相关性**，发现随着肿瘤的进展BTK表达从高到低（之前COX结果中也可以发现它是个保护因素）。对应标题中的"BTK... a prognostic factor"。

- 又分BTK高低组做差异分析，然后用这次得到的差异基因做了GSEA。**发现BTK高组中富集到了一些免疫通路，低组富集到了一些代谢通路，预示着肿瘤微环境在进展过程中发生了一些变化**（对应标题里的"BTK... an indicator for tumor microenvironment remodeling")。

- 最后做CIBERSORT看22种免疫细胞的比例差异（，彩虹图和相关性热图是水图）。又看BTK表达量高低组中这22种免疫细胞的差异，发现了一些有显著差异的细胞类型。再把BTK表达量与22种免疫细胞各自的相对比例每种都做了一个相关性分析，找到显著相关的细胞类型，最后取了交集，发现交集里面取到的这8种免疫细胞就是跟BTK关系比较大的。

## 完结撒花！！！！