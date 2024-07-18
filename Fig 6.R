###
logTPM = read.csv("/share/home/wanglab/wanglab_lh/Bulk_RNAseq/HIVID/logTPM_dereplicated.csv")
rownames(logTPM) = logTPM[,1]
logTPM = logTPM[,-1]

hpv_burden = read.csv("~/CyTOF/HIVID_analysis/HIVID_data/HIVID_burden_samples.csv")
hpv_burden = hpv_burden[,-1]

hpv_burden$group[hpv_burden$group %in% c("HPV_integration_burden_low")] = "Low_Burden"
hpv_burden$group[hpv_burden$group %in% c("HPV_integration_burden_high")] = "High_Burden"

###
data = t(logTPM) %>% as.data.frame()
data$ID = rownames(data)
data
###
data_plot = merge(data,hpv_burden[,1:3],by = "ID")
data_plot = data_plot[,c("ID","TAP1","TAP2","HLA-A","HLA-B","HLA-C",
                         "CXCL9","CXCL10","group")]
data_plot

library(reshape2)
###
plot <- melt(data_plot[,c(2:18)],id = "group")
plot$group = factor(plot$group,levels = c("Low_Burden","High_Burden"))

p = ggboxplot(plot, x = "variable",bxp.errorbar = F,fill = "group",
              y = "value", color = "group",size = 1,add = "none",
              
              xlab = "Genes",ylab = "Expressed level(log2(TPM))",palette = c("#123f81","#cf4728"))  #+stat_compare_means(aex)#+geom_jitter(aes(colour = Cluster))
p = p + stat_compare_means(aes(group= group),method = "wilcox.test",label = "p.format",label.x = 8,size =5) +
  theme(legend.position = 'right',
        axis.text = element_text(size = 20),
        axis.text.x = element_text(size = 20,angle = 30, hjust = 1)) +
  theme(panel.grid.major = element_blank(), #????????
        panel.grid.minor = element_blank(), #????????
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"), #
        panel.background = element_rect(fill = 'white'), #????ɫ
        plot.background= element_rect(fill="white"))

options(repr.plot.width = 24,repr.plot.height = 6)
p

pdf(file = "HLA-I_CXCL_genes_Burden.pdf", width = 20, height = 7)
par(cex = 0.1)
par(mar = c(0,4,2,0))
p
dev.off()

####
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(infercnv)
library(Seurat)
library(copykat)

###
sc = readRDS("~/Analysis/CC_single cell/Tissue/sample_CC_pc30.rds")

tumor = subset(sc,CellType %in% c("Tumor_cell"))
###
objects <- SplitObject(tumor, split.by = "orig.ident")
for (i in 1:length(objects)) {
  objects[[i]] <- NormalizeData(objects[[i]], verbose = FALSE)
  objects[[i]] <- FindVariableFeatures(objects[[i]], selection.method = "vst", 
                                       nfeatures = 2000, verbose = FALSE)
}

sce.anchors <- FindIntegrationAnchors(object.list = objects, dims = 1:10)
objects <- IntegrateData(anchorset = sce.anchors, dims = 1:10)

objects <- ScaleData(objects)
objects <- RunPCA(objects, npcs = 30, verbose = FALSE)

objects <- RunUMAP(objects, reduction = "pca", dims = 1:30)
objects <- FindNeighbors(objects, dims = 1:30)
objects <- FindClusters(objects, resolution = 1)
####
tumor = objects

### Umap plot
library(ggplot2)
my_color = c("#ED593B","#A6CEE3","#79B79A","#CC9B7F","#FBB268",
             "#7C7FC2")

my_color = c("#cf4728","#123f81")

umap = tumor@reductions$umap@cell.embeddings  %>% #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = tumor@meta.data$Integration_Burden) #注释后的label信息 ，改为cell_type


p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  geom_point(size = 1.5 , alpha = 1)  +  
  scale_color_manual(values = my_color)

p2 <- p  +
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))
p3 <- p2 +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=10), #设置legend标签的大小
    legend.key.size=unit(0.5,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5))) #设置legend中 点的大小 

options(repr.plot.width = 7,repr.plot.height = 5)
p3

pdf(file = "~/Analysis/CC_single cell/Tissue/Umap_scRNA_IB.pdf", width = 5.5, height = 4)
par(cex = 0.1)
par(mar = c(0,4,2,0))
p3
dev.off()

### CopyKAT
Fibro = readRDS("~/Analysis/CC_single cell/Tissue/CC_Fibro_qualifued.rds")
Epi = readRDS("~/Analysis/CC_single cell/Tissue/tumor_qualifued.rds")

###
sc = subset(sc,cells = c(rownames(Fibro@meta.data),
                         rownames(Epi@meta.data)),invert = T)

sample_name = as.data.frame(table(sc$Sample))
sample_name
###
for (i in 1:10) {
  test = subset(sc_sub,Sample %in% sample_name[i,1])
  counts = as.matrix(test@assays$RNA@counts)
  
  setwd("/share/home/wanglab/wanglab_lh/Analysis/CC_single cell/Tissue/Copykat")
  
  dir.create(paste0(sample_name[i,1]))
  
  setwd(paste0("/share/home/wanglab/wanglab_lh/Analysis/CC_single cell/Tissue/Copykat/",sample_name[i,1]))
  
  sc_cnv = copykat(rawmat = counts,ngene.chr = 5,sam.name = sample_name[i,1],KS.cut = 0.05, n.cores = 15)
  
  saveRDS(sc_cnv,file = paste0(sample_name[i,1],"_copycat.csv"))
}
