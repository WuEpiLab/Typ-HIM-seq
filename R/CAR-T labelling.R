#----- title: "CAR-T DNA_barcode" -----------------------------------------
#----- prepare R environment ---------------------------------------------
library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(magrittr)
library(pheatmap)
library(tidyverse)
library(cowplot)
library(gridExtra)
library(reticulate)
library(uwot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(Matrix)

#----- set environment ----------------------------------------------------
Sys.setenv(LANGUAGE = "en")
set.seed(123)

#----- loading data ------------------------------------------------------
jia_sc <- Read10X(data.dir = "rawdata/filtered_feature_bc_matrix")
jia_sc = CreateSeuratObject(jia_sc,project = "jia_sc",min_cells = 1,min_ngene = 500,max_ngene = inf)
jia_sc[["percent.mt"]] <- PercentageFeatureSet(jia_sc, pattern = "^MT-")

# vlnplot before control
p1 <- VlnPlot(jia_sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(jia_sc, features =  "nCount_RNA") +
  scale_y_continuous(limits = c(0, 5000),breaks = c(0, 1000, 2000))

ggsave(p1,filename = "result/01.png")

# vlnplot after control
jia_sc <- subset(jia_sc,subset = nCount_RNA > 2000 & percent.mt < 25)
p1 <- VlnPlot(jia_sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "orig.ident")
ggsave(p1,filename = "result/02.png")

#----- seurat Standard pre-processing workflow ----------------------------
jia_sc <- NormalizeData(jia_sc,normalization.method = "LogNormalize", scale.factor = 10000)
jia_sc <- FindVariableFeatures(jia_sc, selection.method = "vst", nfeatures = 2000)

# plot variable features with and without labels
top10 <- head(VariableFeatures(jia_sc), 10)
p1 <- VariableFeaturePlot(jia_sc)
p2 <- LabelPoints(plot = p1, points = top10, repel = TRUE)
plotc = (p1|p2)
ggsave(plotc, filename = "result/03.png",width = 12,height = 6)


all.genes <- rownames(jia_sc)
jia_sc <- ScaleData(jia_sc, features = all.genes)
jia_sc <- RunPCA(jia_sc, verbose = T)

# checking PCA cluster heatmap
# DimHeatmap(jia_sc, dims = 1:20, cells = 500, balanced = TRUE)
p1 <- ElbowPlot(object = jia_sc, ndims = 40)
ggsave(p1, filename = "result/04.png",width = 12,height = 6,bg = "white")
# using ElbowPlot result
jia_sc <- FindNeighbors(jia_sc, reduction = "pca", dims = 1:8)

#Identify suitable resolution
library(clustree)
jia_sc <- FindClusters(object = jia_sc,resolution = c(seq(.1,1.6,.2)))
head(jia_sc@meta.data)

p1 <- clustree(jia_sc@meta.data, prefix = "RNA_snn_res.")
ggsave(p1, filename = "result/05.png",width = 6,height = 12)

#-------resolution decide using the number of cluster
#-------dims decide using the number of PCA
jia_sc <- RunUMAP(jia_sc,dims = 1:8, n.neighbors = 30, min.dist = 0.3)
jia_sc <- RunTSNE(jia_sc,dims = 1:8, perplexity = 30)

p1 <- DimPlot(jia_sc, reduction = "umap",label = T ,group.by = "orig.ident")
p2 <- DimPlot(jia_sc, reduction = "tsne", label = T ,group.by = "orig.ident")
p3 <- DimPlot(jia_sc, reduction = "umap",label = F,group.by = "cell_type") +
  scale_color_manual(values=c('B_cell'='#6bcdff','mix'='#fab838','T_cell'='#9ecb41'))

p4 <- DimPlot(jia_sc, reduction = "tsne", label = T,group.by = "seurat_clusters")

ggsave(p1,filename = "result/06.png", width=8, height=6)
ggsave(p2,filename = "result/07.png", width=8, height=6)
ggsave(p3,filename = "result/08.pdf", width=5, height=4)
ggsave(p4,filename = "result/09.png", width=8, height=6)


saveRDS(jia_sc,file = "result/jia_sc.rds")


jia_sc <- readRDS(file = "result/jia_sc.rds")

#----- dotplot -----------------------------------------------------------------
B_cell <- c("CD79A","MZB1", "MS4A1","CD19", "CD79B")
T_cell <- c("CD2", "CD3D", "CD3E", "CD3G", "CD8A")
cluster_markers <- c("CD79A","MZB1", "MS4A1","CD19", "CD79B","CD2", "CD3D", "CD3E", "CD3G", "CD8A","CD4")
p1 <- DotPlot(jia_sc,features = cluster_markers,group.by = "RNA_snn_res.0.1") + coord_flip()
p1
ggsave(p1,filename = "result/12.png", width = 8, height = 10, bg = "white")





#----- identify cell type ------------------------------------------------------
cluster.ids <- c("0" = "T_cell",
                 "1" = "B_cell",
                 "2" = "mix",
                 "3" = "B_cell")
jia_sc@active.ident <- jia_sc$RNA_snn_res.0.1
jia_sc <- RenameIdents(jia_sc, cluster.ids)
jia_sc$cell_type <- jia_sc@active.ident
p1 <- DimPlot(jia_sc, group.by = "cell_type",label = T)
ggsave(p1,filename = "result/13.png", width = 10, height = 8, bg = "white")

saveRDS(jia_sc,file = "result/jia_sc.rds")


jia_sc <- readRDS(file = "result/jia_sc.rds")




#----- CD19 and DNA_barcode correlation expression  in all cells -------------------
tmp <- jia_sc[['RNA']]@data[c("CD19","DNA_barcode"),]
tmp <- t(tmp)
tmp <- as.data.frame(tmp)
write.csv(tmp,file = "result/CD19_DNA_barcode_allcell.csv")

tmp <- tmp[tmp$CD19 > 0,]
tmp <- tmp[tmp$DNA_barcode > 0,]
genelist <- colnames(tmp)
p1<- ggscatter(tmp, x = "CD19", y = "DNA_barcode",add ="reg.line", conf.int = FALSE,
               cor.coef = TRUE, cor.method ="pearson")

ggsave(p1,filename = "result/14.png", width = 8, height = 8, bg = "white")



table(jia_sc[["RNA"]]@counts[rownames(jia_sc[["RNA"]]@counts) == "CD19",],jia_sc$seurat_clusters)

table(jia_sc[["RNA"]]@counts[rownames(jia_sc[["RNA"]]@counts) == "DNA_barcode",],jia_sc$seurat_clusters)



p1 <- FeatureScatter(subset(jia_sc,CD19 > 0 & DNA_barcode > 0),feature1 = "CD19",feature2 = "DNA_barcode",group.by = "cell_type",pt.size = 2) +
  scale_color_manual(values=c('B_cell'='#6bcdff','mix'='#fab838','T_cell'='#9ecb41'))
ggsave(p1,filename = "result/14-1.pdf", width = 5, height = 4, bg = "white")



#----- loading data ------------------------------------------------------
jia_noscal <- Read10X(data.dir = "rawdata/filtered_feature_bc_matrix")
jia_noscal = CreateSeuratObject(jia_noscal,project = "jia_noscal",min_cells = 1,min_ngene = 500,max_ngene = inf)
jia_noscal[["percent.mt"]] <- PercentageFeatureSet(jia_noscal, pattern = "^MT-")

# after control
jia_noscal <- subset(jia_noscal,subset = nCount_RNA > 2000 & percent.mt < 25)

#----- layout DNA_barcode
DNA_barcode <- jia_noscal[["RNA"]]@counts[rownames(jia_noscal[["RNA"]]@counts) == "DNA_barcode",]
dim(jia_noscal)
tail(rownames(jia_noscal))


#----- add DNA_barcode -------------------------------------------------------------
jia_sc[['RNA']]@counts <- rbind(jia_sc[['RNA']]@counts,DNA_barcode)
jia_sc[['RNA']]@data <- rbind(jia_sc[['RNA']]@data,DNA_barcode)
jia_noscal <- jia_sc
tail(jia_noscal@assays$RNA@data)


saveRDS(jia_noscal,file = "result/jia_noscal.rds")



#----- Bcell apatmer make groups -----------------------------------------------
Bcell_noscal <- jia_noscal[,jia_noscal@meta.data$cell_type == "B_cell"]
two <- colnames(subset(Bcell_noscal,DNA_barcode > 1))
one <- colnames(subset(Bcell_noscal,DNA_barcode < 2 & DNA_barcode > 0))

group <- ifelse(colnames(Bcell_noscal) %in% two,"two",
                ifelse(colnames(Bcell_noscal) %in% one,"one","non"))
table(group)
Bcell_noscal@meta.data$group <- group
head(Bcell_noscal)

saveRDS(Bcell_noscal,file = "result/Bcell_noscal.rds")

#----- DNA_barcode significant gene ------------------------------------------------
two_non <- FindMarkers(Bcell_noscal,
                       ident.1 = "two",
                       ident.2 = "non",
                       group.by = "group",
                       test.use = "wilcox",
                       logfc.threshold = 0)

write.csv(two_non,file = "result/two_non.csv")


# volcano plot
# Set the threshold
# add gene names in excel
two_non <- read.csv(file = "result/two_non.csv")

two_non <- two_non[rownames(two_non) != "DNA_barcode",]


logFC_cutoff <- 0.5
P_Value_cutoff <- 0.05

two_non[which(two_non$avg_log2FC >= logFC_cutoff & two_non$p_val < P_Value_cutoff),'Change'] <- 'up'
two_non[which(two_non$avg_log2FC <= -logFC_cutoff & two_non$p_val < P_Value_cutoff),'Change'] <- 'down'
two_non[which(abs(two_non$avg_log2FC) <= logFC_cutoff | two_non$p_val >= P_Value_cutoff),'Change'] <- 'none'


p1 <- ggplot(data = two_non, aes(x = avg_log2FC, y = -log10(p_val), color = Change)) +
  geom_point(alpha = 0.4, size=2) +  # plot scatter diagram
  scale_color_manual(values = c('red', 'black', 'blue'), limits = c('up', 'none', 'down')) +  # Customize the color of the point
  labs(x = 'avg_log2FC', y = '-log10(p_val)', title = 'two vs non', color = '') +  # axis titles
  theme(plot.title = element_text(hjust = 0.5, size = 14), panel.grid = element_blank(), # background color, grid lines, legend, etc
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 3, color = "black") +  # adding intercept lines
  geom_hline(yintercept = -log10(P_Value_cutoff), lty = 3, color = "black")
#+ xlim(-5, 5) + ylim(0, 15)  # limiting axis scales

library(ggrepel)
two_non$X <- rownames(two_non)
p2 <- p1 + geom_text_repel(data = subset(two_non,abs(p_val) < P_Value_cutoff & abs(avg_log2FC) > logFC_cutoff),
                           aes(label = X),size = 4,col = 'black',box.padding = unit(0.5, "lines"),
                           point.padding = unit(1, "lines"), segment.color = "black", show.legend = FALSE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 100))

ggsave(p1,filename = "result/two_non_1.pdf", width = 8, height = 6)






















