### Cibersort
library(tibble)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library('e1071')             #必须加载，因为后面反卷积的算法是基于这个包里的数据
source("cibersort.R") 

###
data = read.csv("RNAseq.csv")
rownames(data) = data[,1]

write.table(data,"ICB_rnaseq.txt",sep = "\t",col.names = T,row.names = T)   
### run cibersort
results=CIBERSORT("~/GSE_data/Nature_Genetics_CSCC/cibersort_sig.txt","ICB_rnaseq.txt", perm=1000, QN=F)   #还是用那一个test.txt表达矩阵，同样运行

write.csv(results,file = "cibersort_ICB_NG.csv")



####
results = read.csv("cibersort_ICB_NG.csv")
rownames(results) = results[,1]

data_plot = subset(results,Treatment_time %in% c("T1"))
p = ggplot(data_plot, aes(x = Response, y = CD39P_CD8, fill = Response)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2) + 
  scale_fill_manual(values= c("#123f81","#cf4728"))

p = p + stat_compare_means(aes(group= Response),method = "t.test",label = "p.format",label.x = 1.5,size =5) +
  theme(legend.position = 'right',
        axis.text = element_text(size = 20),
        axis.text.x = element_text(size = 20,angle = 30, hjust = 1)) +
  theme(panel.grid.major = element_blank(), #????????
        panel.grid.minor = element_blank(), #????????
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"), #
        panel.background = element_rect(fill = 'white'), #????ɫ
        plot.background= element_rect(fill="white"))
p

pdf(file = "CD39P_pCR vs_non-pCR.pdf",width = 5, height = 7)
par(cex = 1)
par(mar = c(0,4,2,0))
p
dev.off()

###
results$Treatment_time[results$Treatment_time %in% c("T1")] = "Pre"
results$Treatment_time[results$Treatment_time %in% c("T2")] = "Post"

###
non_pCR = subset(results,Response %in% c("non-pCR"))
p = ggpaired(non_pCR, x="Treatment_time", y="CD39P_CD8", fill="Treatment_time",id = "Patient",
             add="jitter",line.color = "gray", line.size = 0.5,
             palette=c("#123f81","#cf4728"),
             xlab=" ", 
             ylab="Percentage（by CIBERSORT）", title = "non-pCR",
             legend.title=" ",show.legend = F) + 
  stat_compare_means(method = "t.test",paired = TRUE, comparisons=list(c("Pre", "Post"))) +#配对t检验
  theme(legend.position = 'none')#去掉legend
p
pdf(file = "CD39P_CD8_non-pCR.pdf",width = 4, height = 7)
par(cex = 1)
par(mar = c(0,4,2,0))
p
dev.off()

###
pCR = subset(results,Response %in% c("pCR"))
p = ggpaired(pCR, x="Treatment_time", y="CD39P_CD8", fill="Treatment_time",id = "Patient",
             add="jitter",line.color = "gray", line.size = 0.5,
             palette=c("#123f81","#cf4728"),
             xlab=" ", 
             ylab="Percentage（by CIBERSORT）", title = "pCR",
             legend.title=" ",show.legend = F) + 
  stat_compare_means(method = "t.test",paired = TRUE, comparisons=list(c("Pre", "Post"))) +#配对t检验
  theme(legend.position = 'none')#去掉legend
p
pdf(file = "CD39P_CD8_pCR.pdf",width = 4, height = 7)
par(cex = 1)
par(mar = c(0,4,2,0))
p
dev.off()


#### pCR score
library(tibble)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggpubr)
library(WGCNA)
library(singscore)
library(tidyverse)
library(GSEABase)

####
info = readxl::read_xlsx("~/41588_2023_1570_MOESM3_ESM.xlsx",11) %>% as.data.frame()
info = info[-1,]
colnames(info) = info[1,]
info = info[-1,]

data = read.csv("~/RNAseq.csv")
rownames(data) = data[,1]
data = data[,100:141]

samples = colnames(data) %>% as.data.frame()
colnames(samples) = "ID"

samples$Patient =str_split(samples$ID,'_',simplify = T)[,1]
samples$stage =str_split(samples$ID,'_',simplify = T)[,2]

samples = merge(samples,info[,c(1,4)],by = "Patient")
samples_Pre = subset(samples,stage %in% c("T1"))

###
expMatrix <- data
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpms <- apply(expMatrix,2,fpkmToTpm)
tpms[1:3,]
colSums(tpms) ### 转换为tpms

###
tpms = tpms[,samples_Pre$ID]

group_list <- samples_Pre$Response
group_list = as.factor(group_list)
group_list = factor(group_list,levels = c("non-pCR","pCR"))

###
#表达矩阵数据校正
exprSet <- tpms
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
#判断数据是否需要转换
exprSet <- log2(exprSet+1)
#差异分析：
dat <- exprSet
design=model.matrix(~factor( group_list ))
fit=lmFit(dat,design) #线性拟合
fit=eBayes(fit) ##Apply empirical Bayes smoothing to the standard errors
options(digits = 4)
topTable(fit,coef=2,adjust='BH')
bp=function(g){
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
deg=topTable(fit,coef=2,adjust='BH',number = Inf)
head(deg) 
#save(deg,file = 'deg.Rdata')

###
logFC_cutoff <- 1
p <- 0.01
k1 = (deg$P.Value< p)&(deg$logFC < -logFC_cutoff)
k2 = (deg$P.Value < p)&(deg$logFC > logFC_cutoff)
deg$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(deg$change)
head(deg)

####
deg_sig = subset(deg,change %in% c("UP","DOWN"))
write.csv(deg_sig,file = "pCR_DEG_sig.csv")


###
#创建GeneSet对象
deg_up = subset(deg,change %in% c("UP"))
deg_down = subset(deg,change %in% c("DOWN"))

up_sig = GeneSet(rownames(deg_up), setName = 'pCR score', geneIdType = SymbolIdentifier())
up_sig
down_sig = GeneSet(rownames(deg_down), setName = 'non-pCR score', geneIdType = SymbolIdentifier())
down_sig

###
#根据基因表达对基因排序
eranks = rankGenes(expr,tiesMethod = "min")
#计算得分
pCR_score = simpleScore(eranks, up_sig,downSet = down_sig,centerScore = FALSE)
head(pCR_score)

write.csv(pCR_score,file = "~/Our_cohort_ICB_signature_singscore.csv")