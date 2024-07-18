###
library(dplyr)
library(tidyverse)
library(patchwork)
library(devtools)

seq <- read.csv("HIVID_Statistic.csv")[,-1]

sample <- as.data.frame(table(seq$Sample))
sample <- sample$Var1

HIVID_HIS_1 <- as.list(1:ncol(seq))
library(stringi)
library(dplyr)
###HIVID_genes_high_percentages
for (i in 1:length(sample)) {
  subset_result <- subset(seq,Sample== sample[i])
  subset_genes <- as.data.frame(table(c(seq$seq1_gene,seq$seq2_gene)))[,1] %>% as.character() 
  subset_genes <- stri_remove_empty(subset_genes)
  HIVID_HIS_1[[i]] <- subset_genes
}

###
dataTemp2 <- do.call(cbind, 
                     lapply(lapply(HIVID_HIS_1, unlist), `length<-`, 
                            max(lengths(HIVID_HIS_1))))
colnames(dataTemp2) <- sample
write.csv(dataTemp2,file = "HIVID_genes in samples.csv")

dataTemp2[is.na(dataTemp2)] = c("")

HIVID_genes <- as.data.frame(table(dataTemp2))

HIVID_genes$Percent <- HIVID_genes[,2]/ncol(seq)
HIVID_genes <- as.data.frame(HIVID_genes)
HIVID_genes <- HIVID_genes[-1,]
write.csv(HIVID_genes,file = "HIVID_gene integration frequency_Statistic.csv")

HIVID_high_genes_seq1 = subset(seq,seq1_gene %in% HIVID_high_genes$dataTemp2)
HIVID_high_genes_seq2 = subset(seq,seq2_gene %in% HIVID_high_genes$dataTemp2)

write.csv(HIVID_high_genes_seq1,file = "Cicros_Plot/HIVID_high_genes_seq1.csv")
write.csv(HIVID_high_genes_seq2,file = "Cicros_Plot/HIVID_high_genes_seq2.csv")
write.csv(HIVID_high_genes,file = "Cicros_Plot/HIVID_high_genes.csv")


####
library(circlize)

circos.info()

cytoband.file = system.file(package = "circlize", "extdata", "cytoBandIdeo.txt")
cytoband.file = system.file("cytoBandIdeo.txt")

cytoband.df = read.table("cytoBandIdeo.txt", colClasses = c("character", "numeric",
                                                            "numeric", "character", "character"), sep = "\t")


cytoband.df = subset(cytoband.df,cytoband.df$V1 %in% c("chr1" , "chr2",  "chr3" , "chr4",  "chr5",  "chr6" , "chr7",  "chr8",
                                                       "chr9"  ,"chr10" ,"chr11","chr12" ,"chr13" ,"chr14", "chr15", "chr16" ,
                                                       "chr17" ,"chr18" ,"chr19" ,"chr20", "chr21" ,"chr22","chrX"))
chromosome <- c("chr1" , "chr2",  "chr3" , "chr4",  "chr5",  "chr6" , "chr7",  "chr8",
                "chr9"  ,"chr10" ,"chr11","chr12" ,"chr13" ,"chr14", "chr15", "chr16",
                "chr17" ,"chr18" ,"chr19" ,"chr20", "chr21" ,"chr22","chrX")
### HIVID info plot
###
circos.clear()
circos.par("start.degree" = 90) ### 摆正位置
circos.initializeWithIdeogram(cytoband.df)
options(repr.plot.width = 30,repr.plot.height = 8)
###
library(rtracklayer)
gtf_data_1 = read.csv("gencode.v42.basic.annotation_rm_duplication.csv")[,-1]
###
colnames(HIVID_genes)[1] = c("gene_name")
HIVID_genes_sites = merge(HIVID_genes,gtf_data_1,by = "gene_name")

###

dir.create("Cicros_Plot")
write.csv(HIVID_genes_sites,file = "Cicros_Plot/HIVID_site_cicrosplot.csv")


a = HIVID_genes_sites1
a$end = as.numeric(a$end)
is.na(a$end) %>% table()
which( is.na (a$end))

for(k in 1:nrow(HIVID_genes_sites)) {
  circos.lines(mean(a[k,5],a[k,6]), -a[k,2],sector.index = a[k,4],
               type = "h",straight = FALSE, col = "#941b0c",baseline = "bottom")
}

####
a = log10(HIVID_results$Sum) %>% as.data.frame()
colnames(a) = "value"
a$variable = "regions"
a$group = "HPV+"
p = ggboxplot(a, x = "variable",bxp.errorbar = F,fill = "value",
              y = "value", color = "group",size = 1,add = "jitter",width = 0,xlab = "all hpv-integration-reads",ylab = "log10(reads)",palette = c("#E71F19","#cf4728"))+coord_flip()#+stat_compare_means(aex)#+geom_jitter(aes(colour = Cluster))
pdf(file = "HIVID_reads_distribution.pdf", width = 10, height = 17)
par(cex = 0.1)
par(mar = c(0,4,2,0))
p
dev.off()

#### Fig2C
HIVID_results = read.csv("HIVID_Statistic.csv")
HIVID_results = as.data.frame(HIVID_results)
HIVID_results = HIVID_results[,-1]

options(repr.plot.width = 9,repr.plot.height = 7)
plot(density(log10(HIVID_results$Sum)), col= "#13956A",ylim = c(0,0.6),lty=1,lwd =3,cex.axis =2,cex.lab=2,cex.main= 2,
     xlab = "Expression level(log10(reads))", main = "HPV integration reads distribution") +
  abline(v = 2.532410)

###Fig2D
plot = read.csv("HPV16_HIVID_statistic.csv")
plot = plot[,-1]

plot$Group = factor(plot$Group,levels = c("Expected","Low_Burden","High_Burden"))
plot$HPV_gene = factor(plot$HPV_gene,levels = c("LCR","E6","E7","E1","E2","E4","E5","L2","L1"))
p1 <- ggplot(plot, aes(HPV_gene, weight = Count, fill = Group)) +
  geom_hline(yintercept = seq(0.1, 0.25, 0.5), color = 'gray') +
  geom_bar(color = "black", width = 0.8, position = 'dodge') +
  labs(x = "HPV breakpoints",y = 'Frequency (%)') +
  scale_fill_manual(values = c("Gray","#123f81","#cf4728")) +
  theme_classic()

pdf(file = "Compared_Freq of HPV break point among expected_Low and high_Burden.pdf", width = 9, height = 5)
par(cex = 0.1)
par(mar = c(0,4,2,0))
p1
dev.off()

####




####
#### MCPcounter data
MCPcounter = read.csv("~/CyTOF/Fig1/MCPcounter_Our cohort.csv")
rownames(MCPcounter) = MCPcounter[,1]
MCPcounter = MCPcounter[,-1]
MCPcounter = t(MCPcounter) %>% as.data.frame()
MCPcounter$ID = rownames(MCPcounter)

MCPcounter = merge(HIVID_Burden,MCPcounter,by = "ID")

###
#data_plot = subset(plot,HPV_State %in% c("Yes"))
library(reshape2)
#plot = subset(plot,HPV_Integration %in% c("Positive"))

data_plot <- melt(MCPcounter[,c(2,4,5,7:13)] ,id= "group")
#data_plot <- subset(data_plot,variable %in%  c("T07","T08","T09","T07"))
data_plot$group <- factor(data_plot$group,levels = c("Low_Burden","High_Burden"))
p = ggboxplot(data_plot, x = "variable",bxp.errorbar = F,fill = "value",
              y = "value", color = "group",size = 1,add = "jitter",xlab = "Cell_Types",ylab = "Expression_level",
              palette = c("#123f81","#cf4728")) +
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)
  ) #+stat_compare_means(aex)#+geom_jitter(aes(colour = Cluster))
p = p + stat_compare_means(aes(group= group),method = "wilcox.test",label = "p.signif",label.x = 8,size =5)
options(repr.plot.width = 20,repr.plot.height = 9)

p+
  theme(
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)
  )

pdf(file = "~/CyTOF/HIVID_analysis/High_Low_Burden_MCPcounter.pdf", width = 12, height = 8)
par(cex = 0.1)
par(mar = c(0,4,2,0))
p
dev.off()
