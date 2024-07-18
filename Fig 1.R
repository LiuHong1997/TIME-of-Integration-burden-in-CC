### Fig 1
#loading packages
library(ggplot2)
library(tidyverse)
library(patchwork)
library(survival)
library(survminer)
library(ggsci)
library(RColorBrewer)
library(ComplexHeatmap)
library(DESeq2)

######Fig1 B
info = readxl::read_xlsx("~/CyTOF/metadata.xlsx",1) %>% as.data.frame()
head(info)

HPV_type_plot = table(info$HPV_types) %>% as.data.frame()
HPV_type_plot = HPV_type_plot[order(HPV_type_plot$Freq,decreasing = T),]

HPV_type_plot$Var1 = factor(HPV_type_plot$Var1,levels = HPV_type_plot$Var1)

####
p1 <- ggplot(HPV_type_plot, aes(Var1, weight = Freq, fill = Var1)) +
  geom_hline(yintercept = seq(0.1, 0.25, 0.5), color = 'gray') +
  geom_bar(color = "black", width = 0.8, position = 'dodge') +
  labs(x = "HPV breakpoints",y = 'Frequency (%)') +
  scale_fill_manual(values = c("#ADD487","#EF7C1B","#369939","#2272A9",
                               "#C6AFD2","#643A90","#AB5728","#F6F29B",
                               "#E57FB0","#D5231E","#00A0E9","#999899",
                               "#F19695","#F5BA6E","#A3C9DD")) +
  theme_classic()

pdf(file = "~/CyTOF/Fig1/Freq of HPV Types.pdf", width = 14, height = 6)
par(cex = 0.1)
par(mar = c(0,4,2,0))
p1
dev.off()

####survial analysis
surv = readxl::read_xlsx("~/CyTOF/surv.xlsx",1) %>% as.data.frame()

###
HPV_surv = merge(surv[,-13],info[,c(3,5)],by = "ID") 
PFS = HPV_surv[,c(7,8,13)] %>% na.omit()
PFS$HPV[!PFS$HPV %in% c("HPV Negative")] = c("HPV Positive")
PFS$PFS_Status = as.numeric(PFS$PFS_Status)

fit=survfit(Surv(PFS_Time,PFS_Status) ~ HPV, data = PFS)

#比较高低表达生存差异
diff=survdiff(Surv(PFS_Time,PFS_Status) ~ HPV,data = PFS)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}

#绘制
surPlot=ggsurvplot(fit, 
                   data=PFS,
                   conf.int=F,
                   pval=pValue,
                   pval.size= 5,
                   legend.labs=c("Negative","Positive"),
                   legend.title= "HPV State",
                   xlab="Time(Month)",
                   break.time.by = 6,
                   risk.table.title="",
                   palette=c("#c54629","#0f3c7a"),
                   risk.table=T,
                   risk.table.height=.25)
print(surPlot)

pdf(file="HPV_State_PFS.pdf",onefile = FALSE,width = 6,height =5)
print(surPlot)
dev.off()

### heatmap
results = read.csv("MCPcounter_Our cohort.csv")
rownames(results) = results[,1]
results = results[,-1]

group <- readxl::read_xlsx("~/CyTOF/metadata.xlsx",1) %>% as.data.frame()
rownames(group) <- group[,1]
group <- group[,-1]
results <- results[,rownames(group)]

library(ggsci)
library(RColorBrewer)
library(ComplexHeatmap)
#myPalette = colorRampPalette(brewer.pal(8,"Set1"))(22)

myPalette = c(brewer.pal(8, "Set1"),brewer.pal(10,"Paired"))

#rdbu <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))

results <- as.matrix(results)

column_ha <- HeatmapAnnotation(
  HPV_integration_Reads = anno_barplot(group$HPV_integration_Reads),
  FIGO_Stage =group$FIGO_Stage_new,
  Pathology = group$Pathology,
  HPV_Clade = group$HPV_Clade,
  HPV_type = group$HPV_types,
  
  HPV_Integration = group$HPV_Integration,
  HPV_State = group$HPV_State,
  
  col = list(
    FIGO_Stage = structure(names = unique(group$FIGO_Stage_new), c(brewer.pal(9,"Set3")[4:7])),
    HPV_Clade = structure(names = unique(group$HPV_Clade), brewer.pal(9,"Set1")[2:9]),
    HPV_type = structure(names = unique(group$HPV_type), c(brewer.pal(12,"Paired")[1:12],brewer.pal(9,"Set1")[8:9])),
    Pathology = structure(names = unique(group$Pathology), c(brewer.pal(9,"Set1")[1:4])),
    HPV_Integration = structure(names = unique(group$HPV_Integration), c(brewer.pal(9,"Pastel1")[1:4])),
    HPV_State = structure(names = unique(group$HPV_State), c(brewer.pal(9,"Set2")[1:2]))
  ),  
  annotation_legend_param = list(
    HPV_Clade = list(
      title = "HPV_Clade",
      at = c("No","A7","A9","A5","A6","A10","A11","Mix Infection"),
      labels = c("No","A7","A9","A5","A6","A10","A11","Mix Infection")
    ),
    HPV_Integration = list(
      title = "HPV_Integration",
      at = c("No","Yes","HPV negative","Undetected"),
      labels = c("No","Yes","HPV negative","Undetected")
    )
  ))

results_scale <- apply(results,1,scale)
results_scale <- t(results_scale)
colnames(results_scale) <- colnames(results)
##
p <- Heatmap(results_scale,top_annotation = column_ha,
             cluster_rows = F,
             cluster_columns = F,
             show_column_names = F)

pdf(file = "ComplexHeatmap.pdf", width = 12, height = 9)
par(cex = 0.1)
par(mar = c(0,4,2,0))
p
dev.off()


####DEG analysis
####DEG_analysis
RNA_count_protein = read.csv("/share/home/wanglab/wanglab_lh/Bulk_RNAseq/RNA_fusion/RNA_count_protein.csv")
RNA_count_protein = RNA_count_protein[,-1]

rownames(RNA_count_protein) = RNA_count_protein[,129]
RNA_count_protein = RNA_count_protein[,-129]
group = read.csv("~/CyTOF/Fig1/cluster_new.csv")
##
###Integration Yes vs No
group2 = subset(group,HPV_types %in% c("HPV non-integration","HPV integration"))
group_list <- group2$HPV_types
group_list = as.factor(group_list)
group_list = factor(group_list,levels = c("HPV non-integration","HPV integration"))
exp = RNA_count_protein[,group2$ID]
colData <- data.frame(row.names =colnames(exp), 
                      condition=group_list)
exp <- round(exp)

dds <- DESeqDataSetFromMatrix(
  countData = exp,
  colData = colData,
  design = ~ condition)
dds <- DESeq(dds)
res <- results(dds,alpha=0.05)
resOrdered <- res[order(res$pvalue),] # ????Pֵ????
DEG_yes_no <- as.data.frame(resOrdered)
head(DEG_yes_no)

logFC_cutoff <- 1
p <- 0.05
k1 = (DEG_yes_no$pvalue< p)&(DEG_yes_no$log2FoldChange < -logFC_cutoff)
k2 = (DEG_yes_no$pvalue < p)&(DEG_yes_no$log2FoldChange > logFC_cutoff)
DEG_yes_no$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG_yes_no$change)
head(DEG_yes_no)

###
library(clusterProfiler)
library(org.Hs.eg.db)

degupgene <- subset(DEG_yes_no,change =="UP")
degdowngene <- subset(DEG_yes_no,change =="DOWN")
degupgene <-  bitr(rownames(degupgene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db )
degdowngene <- bitr(rownames(degdowngene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db )

degupgenego <- enrichGO(gene=degupgene$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        keyType = 'ENTREZID',
                        ont = "BP",
                        pAdjustMethod = 'BH',
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE)

dotplot(degupgenego,showCategory=50, title="HPV Integration associated GOpathways")

dir.create("DEG_analysis")
pdf(file = "DEG_analysis/HPV Integration associated GOpathways.pdf", width = 12, height = 25)
par(cex = 0.6)
par(mar = c(0,4,2,0))
dotplot(degupgenego,showCategory=50, title="HPV Integration associated GOpathways")
dev.off()