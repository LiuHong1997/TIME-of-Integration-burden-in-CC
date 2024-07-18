####
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggsci)
library(RColorBrewer)
library(ComplexHeatmap)
library(survival)
library(survminer)
library(ggpubr)

###
PanelA_tsne = read.csv("tsne.csv")

### TSNE-plot
library(ggplot2)
data1 = PanelA_tsne
allcolour1 = c(
  "#818daf","#8ecabb","#109982","#3a5182",
  "#d94a35","#4db1c9","#7c6047","#ae9a84",
  "#0d6ca9","#d9832b","#23804c","#6c95a8",
  "#fad890","#004383","#8c5a97","#f4ab8f",
  "#aa4643","#698fca"
)
p <- ggplot(data1,aes(x= TSNE1 , y = TSNE2 ,color = cluster)) +  geom_point(size = 1 , alpha =1 )  +  scale_color_manual(values = allcolour1)

p2 <- p  +
  theme(panel.grid.major = element_blank(), #????????
        panel.grid.minor = element_blank(), #????????
        panel.border = element_blank(), #?߿?
        axis.title = element_blank(),  #??????
        axis.text = element_blank(), # ?ı?
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #????ɫ
        plot.background=element_rect(fill="white"))


p3 <- p2 +       
  theme(
    legend.title = element_blank(), #ȥ??legend.title a
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=15), #????legend??ǩ?Ĵ?С
    legend.key.size=unit(0.8,'cm') ) +  # ????legend??ǩ֮???Ĵ?С
  guides(color = guide_legend(override.aes = list(size=5))) #????legend?? ???Ĵ?С 



p4 <- p3 + 
  geom_segment(aes(x = min(data1$TSNE1) , y = min(data1$TSNE2) ,
                   xend = min(data1$TSNE1) +50, yend = min(data1$TSNE2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(data1$TSNE1)  , y = min(data1$TSNE2)  ,
                   xend = min(data1$TSNE1) , yend = min(data1$TSNE2) + 50),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(data1$TSNE1) +25, y = min(data1$TSNE2) -5, label = "tSNE1",
           color="black",size = 6, fontface="bold" ) + 
  annotate("text", x = min(data1$TSNE1) -6, y = min(data1$TSNE2) + 25, label = "tSNE2",
           color="black",size = 6, fontface="bold" ,angle=90) 
p4

pdf(file = "Tissue PanelA_CD3neg cluster.pdf",width = 10, height = 6)
par(cex = 1)
par(mar = c(0,4,2,0))
p3
dev.off()

library(ggplot2)
data1 = PanelA_tsne
colnames(PanelA_tsne)[3] = "cluster"
allcolour1 = c("#dc4932","#5061a3")
p <- ggplot(data1,aes(x= TSNE1 , y = TSNE2 ,color = cluster)) +  geom_point(size = 0.5 , alpha =1 )  +  scale_color_manual(values = allcolour1)

p2 <- p  +
  theme(panel.grid.major = element_blank(), #????????
        panel.grid.minor = element_blank(), #????????
        panel.border = element_blank(), #?߿?
        axis.title = element_blank(),  #??????
        axis.text = element_blank(), # ?ı?
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #????ɫ
        plot.background=element_rect(fill="white"))


p3 <- p2 +         
  theme(
    legend.title = element_blank(), #ȥ??legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=15), #????legend??ǩ?Ĵ?С
    legend.key.size=unit(0.8,'cm') ) +  # ????legend??ǩ֮???Ĵ?С
  guides(color = guide_legend(override.aes = list(size=5))) #????legend?? ???Ĵ?С 

p3

pdf(file = "Tissue PanelA_CD3neg cluster_Tumor_ADJ.pdf",width = 7.5, height = 6)
par(cex = 1)
par(mar = c(0,4,2,0))
p3
dev.off()

###heatmap
heatmap_data = read.csv("~/CyTOF/Fig3/Mat_Mat_Scale_panelACD3-.csv")
rownames(heatmap_data) = heatmap_data[,1]
colnames(heatmap_data)[1] = "Cell_types"

heatmap_data = heatmap_data[c("NK01","NK02", "NK03","NK04","NK05","NK06",
                              "B01","B02","B03","B04","M01","M02","M03","M04","M05","M06","M07","M08"),]

for (i in 3:43) {
  heatmap_data[,i] = scales::rescale(heatmap_data[,i],to = c(0,1))
}

library(ComplexHeatmap)

p = Heatmap(heatmap_data[,-c(1,2,3)], name = "scale",border = "black",
            col = colorRampPalette(c("#f2f4f3","#D9CE8B","#f6bd60","#fb5607","#A62235"))(50),
            cluster_rows = F,cluster_columns = F,rect_gp = gpar(col = "white", lwd = 0.1), 
            left_annotation = rowAnnotation(Cell_types =  heatmap_data$Cell_types,
                                            col = list(Cell_types = structure(names = unique(heatmap_data[,1]), c("#d94a35","#4db1c9","#109982","#3a5182",
                                                                                                                  "#818daf","#8ecabb","#7c6047","#ae9a84",
                                                                                                                  "#0d6ca9","#d9832b","#23804c","#6c95a8",
                                                                                                                  "#fad890","#004383","#8c5a97","#f4ab8f",
                                                                                                                  "#aa4643","#698fca"))),
                                            annotation_legend_param = list(Cell_types = list(
                                              title = "Cell_types",
                                              at = c("NK01","NK02", "NK03","NK04","NK05","NK06",
                                                     "B01","B02", "B03","B04","M01","M02",
                                                     "M03","M04","M05","M06","M07"
                                              ),
                                              labels = c("NK01","NK02", "NK03","NK04","NK05","NK06",
                                                         "B01","B02", "B03","B04","M01","M02",
                                                         "M03","M04","M05","M06","M07"
                                              ))))) +  
  rowAnnotation(Percentage = anno_barplot(heatmap_data$Freq, width = unit(2, "cm")))

pdf(file = "~/CyTOF/Fig3/PanelA_CD45CD3neg_mean_scale_heatmap.pdf", width = 12, height = 7)
par(cex = 0.1)
par(mar = c(0,4,2,0))
p
dev.off()

####
PanelA_CD45 = readxl::read_xlsx("MpanelCD45_T_ADJ_Frequency_boxplot.xlsx") %>% as.data.frame()

#### Tumor ADJ
library(reshape2)
library(ggpubr)
data_plot <- melt(PanelA_CD45[,-2],id= "Group")
data_plot$Group <- factor(data_plot$Group,levels = c("ADJ","Tumor"))
data_plot$variable <- factor(data_plot$variable,levels = c("NK01","NK02", "NK03","NK04","NK05","NK06",
                                                           "B01","B02", "B03","B04","M01","M02",
                                                           "M03","M04","M05","M06","M07","M08"))
p = ggboxplot(data_plot, x = "variable",bxp.errorbar = F,fill = "value",
              y = "value", color = "Group",size = 1,add = "jitter",xlab = "Cell_Types",ylab = "Percentage (% CD45)",palette = c("firebrick","#7768ae"))  #+stat_compare_means(aex)#+geom_jitter(aes(colour = Cluster))
p = p + stat_compare_means(aes(group= Group),method = "wilcox.test",label = "p.signif",label.x = 8,size =5) 
p

pdf(file = "new_version/Tumor_ADJ_CD45_Frequency_compare.pdf", width = 15, height = 7)
par(cex = 0.1)
par(mar = c(0,4,2,0))
p
dev.off()


###Fig4_TSNE_plot
TSNE = read.csv("~/CyTOF/Fig4/tsne.csv")
TSNE = TSNE[,45:47]
TSNE$Cell_types = TSNE$cluster


####
library(ggplot2)
data1 = TSNE
allcolour1 = c("#d94a35","#4db1c9","#109982","#3a5182",
               "#818daf","#8ecabb","#7c6047","#ae9a84",
               "#0d6ca9","#d9832b","#23804c","#6c95a8",
               "#fad890","#004383","#8c5a97","#f4ab8f",
               "#aa4643","#698fca","#959526","#f8c9c9"
)
p <- ggplot(data1,aes(x= TSNE1 , y = TSNE2 ,color = Cell_types)) +  geom_point(size = 0.5 , alpha =1 )  +  scale_color_manual(values = allcolour1)

p2 <- p  +
  theme(panel.grid.major = element_blank(), #????????
        panel.grid.minor = element_blank(), #????????
        panel.border = element_blank(), #?߿?
        axis.title = element_blank(),  #??????
        axis.text = element_blank(), # ?ı?
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #????ɫ
        plot.background=element_rect(fill="white"))


p3 <- p2 +         
  theme(
    legend.title = element_blank(), #ȥ??legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=15), #????legend??ǩ?Ĵ?С
    legend.key.size=unit(0.8,'cm') ) +  # ????legend??ǩ֮???Ĵ?С
  guides(color = guide_legend(override.aes = list(size=5))) #????legend?? ???Ĵ?С 



p4 <- p3 + 
  geom_segment(aes(x = min(data1$TSNE1) , y = min(data1$TSNE2) ,
                   xend = min(data1$TSNE1) +50, yend = min(data1$TSNE2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(data1$TSNE1)  , y = min(data1$TSNE2)  ,
                   xend = min(data1$TSNE1) , yend = min(data1$TSNE2) + 50),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(data1$TSNE1) +25, y = min(data1$TSNE2) -5, label = "tSNE1",
           color="black",size = 6, fontface="bold" ) + 
  annotate("text", x = min(data1$TSNE1) -6, y = min(data1$TSNE2) + 25, label = "tSNE2",
           color="black",size = 6, fontface="bold" ,angle=90) 
p3
pdf(file = "Tissue PanelB_CD3pos cluster.pdf",width = 10.5, height = 7.2)
par(cex = 1)
par(mar = c(0,4,2,0))
p3
dev.off()

####
####
PanelB_CD45_ADJ_Tumor = readxl::read_xlsx("~/CyTOF/Fig4/PanelB_CD3pos_CD45_boxplot.xlsx",2) %>% as.data.frame()
plot = PanelB_CD45_ADJ_Tumor[,c(1,3,5:24)]
colnames(plot)[1] = "Group"

### Tumor ADJ
library(reshape2)
library(ggpubr)
data_plot <- melt(plot[,-2],id= "Group")
data_plot$Group <- factor(data_plot$Group,levels = c("ADJ","Tumor"))
p = ggboxplot(data_plot, x = "variable",bxp.errorbar = F,fill = "value",
              y = "value", color = "Group",size = 1,add = "jitter",xlab = "Cell_Types",ylab = "Percentage (% CD45)",palette = c("firebrick","#7768ae"))  #+stat_compare_means(aex)#+geom_jitter(aes(colour = Cluster))
p = p + stat_compare_means(aes(group= Group),method = "wilcox.test",label = "p.signif",label.x = 8,size =5) 
p

###Diffusion Map
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggsci)
library(RColorBrewer)
library(ComplexHeatmap)
library(survival)
library(survminer)
library(ggpubr)
library(igraph)
library(proxy)
library(diffusionMap)
library(Seurat)
library(destiny)
library(reshape2)
library(Biobase)
library(data.table)

#
TSNE = read.csv("~/CyTOF/Fig4/tsne.csv")
TSNE$Cell_types = TSNE$cluster
TSNE$Cell_types = str_split(TSNE$Cell_types,'_',simplify = T)[,1]
TSNE$X = paste0("No_",TSNE$X)
colnames(TSNE)[1] = "Cell_ID"
###
matrix = TSNE[,1:43] %>% t()
colnames(matrix) = matrix[1,]
matrix = matrix[-1,]

cell_info = TSNE[,c(1,44:58)]

a = cell_info$Cell_types %>% as.data.frame()
rownames(a) = rownames(cell_info)
colnames(a) = "color"
a$cell_types = a$color
##
a$color[a$color %in% c('T06')] = "#8ecabb"
a$color[a$color %in% c('T07')] = "#7c6047"
a$color[a$color %in% c('T08')] = "#ae9a84"
a$color[a$color %in% c('T09')] = "#0d6ca9"
a$color[a$color %in% c('T13')] = "#fad890"
a$color[a$color %in% c('T14')] = "#004383"
a$color[a$color %in% c('T18')] = "#698fca"
a$color[a$color %in% c('T19')] = "#959526"

###
###
matrix = TSNE[,1:43] %>% t()
colnames(matrix) = matrix[1,]
matrix = matrix[-1,]

CyTOF = CreateSeuratObject(counts = matrix,meta.data = cell_info)
library(Biobase)
ct <-GetAssayData(object = CyTOF)

ct <- as.ExpressionSet(as.data.frame(t(ct)))
ct
#. Annotations can be accessed directly via ct$column and ct[['column']]. 
ct$celltype <- cell_info$Cell_types
dm <- DiffusionMap(ct,k = 10)

saveRDS(dm,file = "Diffusion_map_dm.rds")

palette(cube_helix(20)) # configure color palette
plot(dm, pch = 20) # pch for prettier points

####

