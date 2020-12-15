library(dplyr)
library(Seurat)
library(Matrix)
library(cowplot)
library(ggplot2)

# Specify data path
data_path <- "dataset/"


# load datsets
countdata <- read.csv(paste0(data_path, "GSE89567_normalized.csv.gz"), header = TRUE, 
                       row.names = 1, as.is = TRUE, stringsAsFactors = FALSE)
metadata <- read.csv(paste0(data_path, "metadata.csv"), header = TRUE, 
                     row.names = 1, as.is = TRUE)

countdata <- Matrix(as.matrix(countdata), sparse = TRUE)

Brain <- CreateSeuratObject(counts = countdata, project ="Brain", 
                            min.cells = 5, min.features = 100, meta.data = metadata)

Brain <- FindVariableFeatures(Brain, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Brain), 10)
top10

plot1 <- VariableFeaturePlot(Brain)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(Brain)
Brain <- ScaleData(Brain, features = all.genes)

# Clustering the cell using PCA & UMAP
Brain <- RunPCA(Brain, features = VariableFeatures(object = Brain))
ElbowPlot(Brain)

Brain <- FindNeighbors(Brain, dims = 1:15)
Brain <- FindClusters(Brain, resolution = 0.5)
DimPlot(Brain, reduction = "pca")

Brain <- RunUMAP(Brain, dims = 1:15)
DimPlot(Brain, reduction = "umap")

# Changing the identity with metadata
Idents(Brain) <- Brain@meta.data$leiden
head(Idents(Brain))
p1 <- DimPlot(Brain, reduction = "umap")

Idents(Brain) <- Brain@meta.data$Major_type
head(Idents(Brain))
p2 <- DimPlot(Brain, reduction = "umap")

plot_grid(p1, p2)

Tumor.markers <- FindMarkers(Brain, ident.1 = "Tumor", min.pct = 0.25)
Macrophage.markers <- FindMarkers(Brain, ident.1 = "Macrophage", min.pct = 0.25)

head(Tumor.markers, n = 10)
head(Macrophage.markers, n = 10)

write.table(Tumor.markers, file="Tumor.markers.txt", sep='\t')
write.table(Macrophage.markers, file="Macrophage.markers.txt", sep='\t')

FeaturePlot(Brain, features = c("'EGFR'", "'CD14'"))

