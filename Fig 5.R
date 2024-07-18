###
library(ggplot2)
library(tidyverse)
library(patchwork)
library(Seurat)
library(Matrix)
library(forcats)
library(reshape2)
library(ggpubr)


###
fs=list.files(pattern = '.h5')
fs
library("hdf5r")

sceList =lapply(fs, function(x){
  # x=fs[1]
  print(x)
  a=Read10X_h5( x )
  a[1:4,1:4] 
  library(stringr)
  (p=str_split(x,'_',simplify = T)[,1])
  sce <- CreateSeuratObject( a ,project = p )
  sce
})

####
sc <- RenameCells(sc, add.cell.id = sc$orig.ident)  ##
mt.genes <- rownames(sc)[grep("^MT-",rownames(sc), ignore.case = TRUE)]
sc$percent.mt <- PercentageFeatureSet(sc, features = mt.genes )
## ͳ?ƺ???????????��??ͬ????Ҫ֪?��??????????????ֹ???"RPS"??"RPL"
rb.genes <- rownames(sc)[grep("^RP[SL]",rownames(sc), ignore.case = TRUE)]
sc$percent.rb <- PercentageFeatureSet(sc, features = rb.genes )

###
sc <- subset(sc, subset = nCount_RNA > 300 & nCount_RNA < 30000 & nFeature_RNA > 600 & nFeature_RNA < 6000 & percent.mt < 15)

###
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb")
group = "orig.ident"
theme.set2 = theme(axis.title.x=element_blank())

plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(CD3_PBMC, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}

violin <- wrap_plots(plots = plots, nrow=2)  
violin

####
library(DoubletFinder)
library(dplyr)

# Define doublet formation rate
doublet_formation_rate <- 0.075
run_doubletfinder <- function(s) {
  # Process normally
  s <- s %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
  
  # Run TSNE clustering
  s <- RunPCA(s, features = VariableFeatures(object = s))
  s <- FindNeighbors(s, dims = 1:12)
  s <- FindClusters(s, resolution = 0.3)
  s <- RunUMAP(s, dims = 1:12)
  
  # pK Identification (no ground-truth) 
  sweep.res.list <- paramSweep_v3(s, PCs = 1:12, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  max_index <- which.max(bcmvn$BCmetric)
  optimal_pK <- as.numeric(as.character(bcmvn[max_index, "pK"]))
  
  # Homotypic Doublet Proportion Estimate
  annotations <- s@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_formation_rate*nrow(s@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  print(s)
  
  # Run DoubletFinder
  s <- doubletFinder_v3(s, PCs = 1:12, pN = 0.25, pK = optimal_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  # Rename column name for consistency
  colnames(s@meta.data)[grep("DF.classifications*", colnames(s@meta.data)) ] <- "DF.classifications"
  print(head(s@meta.data))
  return(s)
}
sc =  SplitObject(sc, split.by = "orig.ident")
seurat_object_list_1 <- run_doubletfinder(sc[[1]])
seurat_object_list_2 <- run_doubletfinder(sc[[2]])
seurat_object_list_3 <- run_doubletfinder(sc[[3]])

sc <- merge(x = seurat_object_list_1,y = c(
  seurat_object_list_2,
  seurat_object_list_3
))

sc <- SplitObject(sc, split.by = "orig.ident")

for (i in 1:length(sc)) {
  sc[[i]] <- NormalizeData(sc[[i]], verbose = FALSE)
  sc[[i]] <- FindVariableFeatures(sc[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
}

sce.anchors <- FindIntegrationAnchors(object.list = sc, dims = 1:30)
sc <- IntegrateData(anchorset = sce.anchors, dims = 1:30)

sc <- ScaleData(sc)
sc <- RunPCA(sc, npcs = 30, verbose = FALSE)

sc <- RunUMAP(sc, reduction = "pca", dims = 1:30)
sc <- FindNeighbors(sc, dims = 1:30)
sc <- FindClusters(sc, resolution = 1)


options(repr.plot.width = 12,repr.plot.height = 12)
DimPlot(CD3_PBMC,label = T)


#### TCR clonality compare

###
meta_CD39pos = subset(sc@meta.data,Group %in% c("CD39 positive"))
meta_CD39pos = subset(meta_CD39pos,!Clonotype %in% c("unpaired TCR sequence"))
Cellratio <- prop.table(table(meta_CD39pos$Clonotype,meta_CD39pos$Sample_CAT), margin = 2) %>% as.data.frame()
colnames(Cellratio)[1] = "group"
Cellratio[,1] = as.character(Cellratio[,1])
Cellratio[,2] = as.character(Cellratio[,2])
Cellratio_1 = Cellratio
Cellratio_1$clusters = "CD39_positive"
###
###
meta_CD39neg = subset(sc@meta.data,Group %in% c("CD39 negative"))
meta_CD39neg = subset(meta_CD39neg,!Clonotype %in% c("unpaired TCR sequence"))
Cellratio <- prop.table(table(meta_CD39neg$Clonotype,meta_CD39neg$Sample_CAT), margin = 2) %>% as.data.frame()
colnames(Cellratio)[1] = "group"
Cellratio[,1] = as.character(Cellratio[,1])
Cellratio[,2] = as.character(Cellratio[,2])
Cellratio_2 = Cellratio
Cellratio_2$clusters = "CD39_negative"

###
Cellratio = rbind(Cellratio_1,Cellratio_2)

Cellratio$clusters <- factor(Cellratio$clusters,levels = c("CD39_negative","CD39_positive"))
Cellratio$group <- factor(Cellratio$group,levels = c("n = 1","1 < n <= 5","n > 5"))
p = ggboxplot(Cellratio, x = "group",bxp.errorbar = F,fill = "Freq",
              y = "Freq", color = "clusters",size = 1,add = "jitter",xlab = "Clonotypes",ylab = "Percentage",
              palette = c("#123f81","#cf4728")) 
p = p + stat_compare_means(aes(group= clusters),method = "t.test",label = "p.signif",label.x = 8,size =5) 
p

pdf(file = "~/scRNAseq/CCA-TIL/CD39+_-_TCRclone_compare_all clonotypes.pdf", width = 4, height = 4)
par(cex = 0.1)
par(mar = c(0,4,2,0))
p
dev.off()

####
