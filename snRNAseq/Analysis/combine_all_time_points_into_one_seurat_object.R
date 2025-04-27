# Combine all time points into one seurat object
# Silje Wilhelmsen

# Load necessary libraries
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

# Load individual Seurat objects
mouse_vcm_6h <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_6h_vcm.Rds")
mouse_vcm_12h <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_12h_vcm.Rds")
mouse_vcm_1d <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_1d_vcm.Rds")
mouse_vcm_3d <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_vcm_3d.Rds")
mouse_vcm_1w <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_vcm_1w.Rds")
mouse_vcm_3w <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_vcm_3w.Rds")

# # If already in the environment
# mouse_vcm_6h <- readRDS("C:/Users/Labuser/master_project/snRNAseq/tmp/mouse_6h_vcm.Rds")
# mouse_vcm_12h <- readRDS("C:/Users/Labuser/master_project/snRNAseq/tmp/mouse_12h_vcm.Rds")
# mouse_vcm_1d <- readRDS("C:/Users/Labuser/master_project/snRNAseq/tmp/mouse_1d_vcm.Rds")
# mouse_vcm_3d <- readRDS("C:/Users/Labuser/master_project/snRNAseq/tmp/mouse_vcm_3d.Rds")
# mouse_vcm_1w <- readRDS("C:/Users/Labuser/master_project/snRNAseq/tmp/mouse_vcm_1w.Rds")
# mouse_vcm_3w <- readRDS("C:/Users/Labuser/master_project/snRNAseq/tmp/mouse_vcm_3w.Rds")


# Merge Seurat objects
mouse_vcm_all_time_points <- merge(mouse_vcm_6h, 
                                   y = list(mouse_vcm_12h, 
                                            mouse_vcm_1d, 
                                            mouse_vcm_3d, 
                                            mouse_vcm_1w, 
                                            mouse_vcm_3w), 
                                   project = "mouse_full_merge_project", 
                                   merge.data = TRUE)
mouse_vcm_all_time_points <- NormalizeData(mouse_vcm_all_time_points, normalization.method = "LogNormalize", scale.factor = 10000)
mouse_vcm_all_time_points <- FindVariableFeatures(mouse_vcm_all_time_points)
mouse_vcm_all_time_points <- ScaleData(mouse_vcm_all_time_points)
mouse_vcm_all_time_points <- RunPCA(mouse_vcm_all_time_points, npcs = 20)
mouse_vcm_all_time_points <- RunUMAP(mouse_vcm_all_time_points, dims = 1:20)


# Save Seurat object in h5Seurat format
SaveH5Seurat(mouse_vcm_all_time_points, filename = "C:/Users/Labuser/snRNAseq/tmp/mouse_all_time_points.h5seurat", overwrite = TRUE)

# Save Seurat object in .Rds format
saveRDS(mouse_vcm_all_time_points, "C:/Users/Labuser/snRNAseq/tmp/mouse_vcm_all_time_points.Rds")



mouse_vcm_all_time_points <- IntegrateLayers(object = mouse_vcm_all_time_points, 
                                              method = RPCAIntegration, 
                                              orig.reduction = "pca", 
                                              new.reduction = "integrated.rpca",
                                              verbose = FALSE)
mouse_vcm_all_time_points_join_layers <- JoinLayers(mouse_vcm_all_time_points)
                                             
# Save Seurat object in h5Seurat format
SaveH5Seurat(mouse_vcm_all_time_points_join_layers, filename = "C:/Users/Labuser/snRNAseq/tmp/mouse_all_time_points_join_layers.h5seurat", overwrite = TRUE)

# Save Seurat object in .Rds format
saveRDS(mouse_vcm_all_time_points_join_layers, "C:/Users/Labuser/snRNAseq/tmp/mouse_vcm_all_time_points_join_layers.Rds")



# Create UMAP plot
dimplot_umap <- DimPlot(mouse_vcm_all_time_points, reduction = "umap", label = TRUE) + 
  ggtitle("UMAP Plot") +
  theme_minimal()  

print(dimplot_umap)

# If you intended to integrate using RPCA, follow this approach:
# 1. Split object by original identity (time points in this example)
mouse_vcm.list <- SplitObject(mouse_vcm_all_time_points, split.by = "orig.ident")

# 2. Normalize and identify variable features for each dataset independently
mouse_vcm.list <- lapply(mouse_vcm.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
  return(x)
})

# 3. Select features for integration
features <- SelectIntegrationFeatures(object.list = mouse_vcm.list)
mouse_vcm.list <- lapply(X = mouse_vcm.list, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
  return(x)
})

# 4. Find integration anchors
anchors <- FindIntegrationAnchors(object.list = mouse_vcm.list, reduction = "rpca", dims = 1:20)

# 5. Integrate data
mouse_vcm_all_time_points <- IntegrateData(anchorset = anchors, dims = 1:20)

# Optional: Run post-integration analysis
# Scale the data again post-integration
mouse_vcm_all_time_points <- ScaleData(mouse_vcm_all_time_points)

# Run PCA again post-integration
mouse_vcm_all_time_points <- RunPCA(mouse_vcm_all_time_points, npcs = 20)

# Run UMAP again post-integration
mouse_vcm_all_time_points <- RunUMAP(mouse_vcm_all_time_points, dims = 1:20)

# Save Seurat object in h5Seurat format
SaveH5Seurat(mouse_vcm_all_time_points, filename = "mouse_all_time_points.h5seurat", overwrite = TRUE)

# Save Seurat object in .Rds format
saveRDS(mouse_vcm_all_time_points, "C:/Users/siljeew/snRNAseq/tmp/mouse_vcm_all_time_points.Rds")

# Create UMAP plot
dimplot_umap <- DimPlot(mouse_vcm_all_time_points, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE) + 
  ggtitle("UMAP Plot") +
  theme_minimal()

# Explanation of Steps:
# 1. Loading Data: Load each of your time-point-specific Seurat objects from .RDS files.
# 2. Merging Data: Merge individual Seurat objects using merge().
# 3. Preprocessing:
#   Normalize the data.
#   Identify variable features.
#   Scale the data.
#   Run PCA.
#   Run UMAP.
# 4. Integration:
#   Split data by original identity (time points).
#   Normalize and identify variable features for each subset independently.
#   Select integration features.
#   Scale data and run PCA for each subset.
#   Find integration anchors using RPCA reduction.
#   Integrate data using the anchors.
#   Optionally, re-scale, re-run PCA, and re-run UMAP for integrated data.
# 5. Saving Data: Save the final integrated Seurat object both as .h5Seurat and .RDS files.


# # This is copied from the individual preprocessing scripts
# ################################################################################
# # Combine all time points into one seurat object
# # Silje Wilhelmsen
# 
# 
# library(Seurat)
# library(SeuratData)
# library(ggplot2)
# library(patchwork)


options(future.globals.maxSize = 48 * 1024^3)  # Setting max size to 48 GiB in case of error due to size.

mouse_6h <- readRDS("C:/Users/Labuser/snRNAseq/tmp/mouse_6h.Rds")
mouse_12h <- readRDS("C:/Users/Labuser/snRNAseq/tmp/mouse_12h.Rds")
mouse_1d <- readRDS("C:/Users/Labuser/snRNAseq/tmp/mouse_1d.Rds")
mouse_3d <- readRDS("C:/Users/Labuser/snRNAseq/tmp/mouse_3d.Rds")
mouse_1w <- readRDS("C:/Users/Labuser/snRNAseq/tmp/mouse_1w.Rds")
mouse_3w <- readRDS("C:/Users/Labuser/snRNAseq/tmp/mouse_3w.Rds")


mouse_all_time_points <- merge(mouse_6h,
                                y = c(mouse_12h,
                                     mouse_1d,
                                     mouse_3d,
                                     mouse_1w,
                                     mouse_3w),
                                project = "mouse_full_merge_project",
                                merge.data = TRUE)

mouse_all_time_points <- NormalizeData(mouse_all_time_points,
                                           normalization.method = "LogNormalize",
                                           scale.factor = 10000)
mouse_all_time_points <- FindVariableFeatures(mouse_all_time_points)
mouse_all_time_points <- ScaleData(mouse_all_time_points)
mouse_all_time_points <- RunPCA(mouse_all_time_points, npcs = 20)
mouse_all_time_points <- RunUMAP(mouse_all_time_points, dims = 1:20)
mouse_all_time_points <- IntegrateLayers(object = mouse_all_time_points,
                                             method = RPCAIntegration,
                                             orig.reduction = "pca",
                                             new.reduction = "integrated.rpca",
                                             verbose = FALSE)


mouse_all_time_points_join_layers <- JoinLayers(mouse_all_time_points)

SaveH5Seurat(mouse_all_time_points_join_layers,
             filename = "mouse_all_time_points_join_layers.h5seurat",
             overwrite = TRUE)

saveRDS(mouse_all_time_points_join_layers, file = "C:/Users/Labuser/snRNAseq/tmp/mouse_all_time_points_join_layers.Rds")
