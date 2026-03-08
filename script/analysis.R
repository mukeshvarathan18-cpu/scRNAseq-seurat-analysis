library(Seurat)

# load dataset (10X format)
data <- Read10X("C:/Users/HP/Documents/scRNAseq")

# check dimension of the matrix (genes, cells)
dim(data)

# create seurat object from count matrix
seurat_1 <- CreateSeuratObject(counts = data)

# to view information about seurat object
seurat_1
#----------------------
### QUALITY CONTROL
#----------------------
# Count the mitochondrial gene expression percentage
seurat_1[["percent.mt"]] <- PercentageFeatureSet(seurat_1, pattern = "^MT-")

# violin plot which represents the QC matrices
VlnPlot(seurat_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

# Filter out the low quality cells
seurat_1 <- subset(
  seurat_1,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 5
)

#-----------------------
### DATA NORMALIZATION
#-----------------------
# Normalize gene expression value
seurat_1 <- NormalizeData(seurat_1)

#--------------------------------------
### IDENTIFYING HIGHLY VARIABLE GENES
#--------------------------------------
# Among the all genes detecting higly variable genes across the cells
seurat_1 <- FindVariableFeatures(seurat_1)

# To view the plot of variable genes
VariableFeaturePlot(seurat_1)

# list out top variable genes
head(VariableFeatures(seurat_1))

#-----------------
### DATA scaling
#-----------------
# to remove the unwanted variation in data
seurat_1 <- ScaleData(seurat_1)

#----------------------------------
### PRINCIPAL COMPONENT ANALYSIS
#----------------------------------
# To perform the PCA
seurat_1 <- RunPCA(seurat_1)

# Number of PCs to use
ElbowPlot(seurat_1)
print(seurat_1[["pca"]], dims = 1:5, nfeatures = 5)
ElbowPlot(seurat_1)

#---------------------
### CLUSTERING CELLS
#---------------------
# To detect the nearest neighbor graph in cells
seurat_1 <- FindNeighbors(seurat_1, dims = 1:10)

# Identify the cluster cells
seurat_1 <- FindClusters(seurat_1, resolution = 0.5)

#-----------------------
### UMAP VISUALIZATION 
#-----------------------
# To run UMAP visualization
seurat_1 <- RunUMAP(seurat_1, dims = 1:10)

# plot clusters on UMAP
DimPlot(seurat_1, reduction = "umap", label = TRUE)

#-------------------------------
### MARKER GENE IDENTIFICATION
#-------------------------------
# Identify the marker gene for each clusters
marker <- FindAllMarkers(seurat_1)

# Top marker genes
head(marker)

# gene expression visualization
FeaturePlot(seurat_1, features = "CD3D")

# top 10 marker gene of each clusters
library(dplyr)
top10 <- marker %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC)

# Heatmap of top marker genes
DoHeatmap(seurat_1, features = top10$gene) + NoLegend()

#----------------------------
### ADDITIONAL VISUALIZATION
#-----------------------------
head(rownames(seurat_1))
FeaturePlot(seurat_1, features = rownames(seurat_1)[1])
DotPlot(seurat_1, features = unique(top10$gene)) + RotatedAxis()

--------------------
### RENAME CLUSTERS
--------------------
# Assigning new cluster names
levels(seurat_1)

# Map name to cluster levels
new.cluster.ids <- c("cluster_1","cluster_2","cluster_3","cluster_4","cluster_5")
names(new.cluster.ids) <- levels(seurat_1)

# Rename cluster identities
seurat_1 <- RenameIdents(seurat_1, new.cluster.ids)

# To visualize the renamed clusters
DimPlot(seurat_1, reduction = "umap", label = TRUE)

# ----------------------
# SAVE RESULTS
# ----------------------

# Save Seurat object
saveRDS(seurat_1, file = "scRNA_project.rds")

# Export marker gene tables
write.csv(marker, "cluster_markers.csv")
write.csv(top10, "top10markers.csv")
