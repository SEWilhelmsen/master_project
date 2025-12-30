# Run preprocessing of spatial transcriptomics data
# https://satijalab.org/seurat/articles/visiumhd_commands_intro#visium-hd-support-in-seurat

# Load libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggplot2)

# Set global max size
options(future.globals.maxSize = 5.0 * 1e9) # Allow 5 GB
load("C:/Users/Labuser/master_project/SCreference_l2_ALLNUCLEI_10112025.RData")


# Load data
### 1 day ### 
### ORAB ###
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D2-TP39-4-1/outs"
# TP394 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# object <- TP394
# 
# # ## 1 day ### 
# # ### SHAM ###
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A2-TP38-4-2/outs"
# TP384 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# object <- TP384
#  
# ## 3 days ### 
# ### ORAB ###
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D1-TP36-1-1/outs"
# TP361 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# object <- TP361
# 
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D1-TP37-4-1/outs"
# TP374 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# 
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D2-TP37-2/outs"
# TP372 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# 
# ## 3 days ### 
# ### SHAM ###
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A1-TP37-5-3/outs"
# TP375  <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# 
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A1-TP36-3/outs"
# TP463 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# 
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A2-TP36-4/outs"
# TP364 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# 
# ## 1 Week ### 
# ### ORAB ###
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D2-TP52-4/outs"
# TP524 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# 
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A2-TP52-3/outs"
# TP523 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# 
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D2-TP46-4/outs"
# TP464 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# 
# ## 1 Week ### 
# ### SHAM ###
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A2-TP46-2/outs"
# TP462 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# 
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D1-TP42-1/outs"
# TP421 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# 
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D1-TP42-2/outs"
# TP422 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# 
# ## 3 Week ### 
# ### ORAB ###
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D1-TP32-3-3/outs"
# TP323 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# 
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0901/GCF0901/A1/outs"
# A1 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# 
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0901/GCF0901/D2/outs"
# D2 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# 
# ## 3 Weeks ### 
# ### SHAM ###
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0901/GCF0901/D1/outs"
# D1 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# 
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A1-TP48-1-1/outs"
# TP481 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
# 
# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A1-TP31-1/outs"
# localdir <- "C:/Users/Labuser/master_project/Spatials/tmp/A1-TP31-1/outs"
# TP311 <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))


# object <- merge(A1, y = c(D1, D2, TP323, TP361, TP364, TP372, TP374, TP375, TP384, TP394, TP421, TP422, TP464, TP481, TP523, TP524))
object <- merge(TP361, y = c(TP372, TP374, TP375, TP384, TP394))

head(object@meta.data$orig.ident)

# can use Spatial.016um to analyze 16um binning as well
DefaultAssay(object) <- "Spatial.008um"
object <- NormalizeData(object)
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, reduction.name = "pca.008um") # Time: 4+ hr 
object <- FindNeighbors(object, reduction = "pca.008um", dims = 1:20, graph.name = "SNN.008um", k.param = 30) # Time: 20 min. dims = 1:30 produced error.
# Reduce k.param for more a more detailed graph.

names(object@graphs)

object <- FindClusters(object, graph.name = "SNN.008um", resolution = 0.6, cluster.name = "seurat_cluster.008um") # Time: 8 hrs
#object <- FindClusters(object, graph.name = "SNN.008um", resolution = 0.005, cluster.name = "seurat_cluster.008um") # Time: 8 hrs. Clusters: 56
# object <- FindClusters(object, graph.name = "SNN.008um", resolution = 0.002, cluster.name = "seurat_cluster.008um") # Time: 8 hrs

object <- RunUMAP(object, reduction = "pca.008um", reduction.name = "umap.008um", dims = 1:20) # Time: 4 hrs


saveRDS(object, file = "C:/Users/Labuser/master_project/snRNAseq/tmp/spatial_object_with_umap_tp361.rds")
object <- readRDS("C:/Users/Labuser/master_project/snRNAseq/tmp/spatial_object_with_umap_tp361.rds")

# Save plot as spatial_object_with_umap_tp...tiff
dim.plot <- DimPlot(object, reduction = "umap.008um", group.by = "seurat_cluster.008um", label = TRUE,
                   repel = T) + NoLegend()
dim.plot

cluster.plot <- SpatialDimPlot(object, group.by = "seurat_cluster.008um", label = FALSE, pt.size.factor = 1.2) +
 theme(legend.position = "right")

cluster.plot




# Make p5 and p6 for 2-3 clusters from the violin plot, and save individual images if possible instead of displaying all in one plot
# can highlight by identity for more detailed visualization
p0 <- SpatialDimPlot(object, cells.highlight = WhichCells(object, expression = seurat_cluster.008um ==
                                                            0), pt.size.factor = 1.2) + NoLegend() + ggtitle("Cluster 0")
p0

p1 <- SpatialDimPlot(object, cells.highlight = WhichCells(object, expression = seurat_cluster.008um ==
                                                            1), pt.size.factor = 1.2) + NoLegend() + ggtitle("Cluster 1")
p1
p2 <- SpatialDimPlot(object, cells.highlight = WhichCells(object, expression = seurat_cluster.008um ==
                                                            2), pt.size.factor = 1.2) + NoLegend() + ggtitle("Cluster 2")
p2
p3 <- SpatialDimPlot(object, cells.highlight = WhichCells(object, expression = seurat_cluster.008um ==
                                                            3), pt.size.factor = 1.2) + NoLegend() + ggtitle("Cluster 3")
p3
p4 <- SpatialDimPlot(object, cells.highlight = WhichCells(object, expression = seurat_cluster.008um ==
                                                            5), pt.size.factor = 1.2) + NoLegend() + ggtitle("Cluster 5")
p4


# Create violin plot of clusters and samples: 
open("C:/Users/Labuser/master_project/Spatials/run_create_violin_plot_of_clusters_and_samples.R")


# Prepare for further visualization
DefaultAssay(object) <- "Spatial.016um"
counts_hd <- object[["Spatial.016um"]]$counts
object_cells_hd <- colnames(object[["Spatial.016um"]])
coords <- GetTissueCoordinates(object)[object_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

rownames(SCreference@counts) <- toupper(rownames(SCreference@counts))
rownames(query@counts) <- toupper(rownames(query@counts))

RCTD <- create.RCTD(query, SCreference, max_cores = 6)
RCTD <- run.RCTD(RCTD, doublet_mode = "full") # this command takes ~15 mins to run

barcodes <- colnames(RCTD@spatialRNA@counts) # list of spatial barcodes
weights <- RCTD@results$weights # Weights for each cell type per barcode

norm_weights <- normalize_weights(weights)
cell_type_names <- colnames(norm_weights) # List of cell types
subset_df <- as.data.frame(t(as.data.frame(norm_weights[1:2,])))
subset_df$celltypes <- rownames(subset_df); rownames(subset_df) <- NULL

i <- 19
plot_puck_continuous(RCTD@spatialRNA, barcodes, norm_weights[,cell_type_names[i]], title =cell_type_names[i], size=1)
# 600 dpi, tiff, 
# Bruk celltype 19 (ventrikulære cm) for å sjekke om det gir mening med fordeling av celletype

plot_dir <- "C:/Users/Labuser/master_project/Spatials/Plots/1 Week/SHAM/"

# Lagre plottet som en TIFF-fil
tiff(paste0(plot_dir, prefix, "plot_puck_", i, ".tiff"), width = 10, height = 10, units = "in", res = 600)
plot_puck_continuous(RCTD@spatialRNA, barcodes, norm_weights[, cell_type_names[i]], title = cell_type_names[i], size = 1)
dev.off()



