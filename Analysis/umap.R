### Create UMAP-plot
# Input:
# Output: visualization of feature expression
# Help: https://satijalab.org/seurat/articles/visualization_vignette 
## Silje Wilhelmsen

# Do UMAP on the all cells and all times -> cardiomyocytes should always be clustered as a cardiomyocyte
# Do UMAP on subtypes of cells, e.g. are there cluster possibilities within the cells?

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)

dimplot_umap <- DimPlot(mouse_vcm_all_time_points, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE) + 
  ggtitle("UMAP Plot") +
  theme_minimal()

# # Display the plot
# print(dimplot_umap_mouse_3w)
# 
# library(Seurat)
# library(ggplot2)
# library(ComplexHeatmap)
# 
# # Assuming 'mouse_3w' is already processed and contains UMAP embeddings.
# # Normalize, find variable features, scale, run PCA, and run UMAP
# mouse_3w <- NormalizeData(mouse_3w, normalization.method = "LogNormalize", scale.factor = 10000)
# mouse_3w <- FindVariableFeatures(mouse_3w)
# mouse_3w <- ScaleData(mouse_3w)
# mouse_3w <- RunPCA(mouse_3w, npcs = 20)
# mouse_3w <- RunUMAP(mouse_3w, dims = 1:20)

# Define your custom colors
custom_colors <- c("ORAB" = "coral", "SHAM" = "black")

dimplot_umap_mouse_3w <- DimPlot(mouse_3w, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE) + 
  ggtitle("UMAP Plot of mouse_3w") +
  theme_minimal()

# Display the plot
print(dimplot_umap_mouse_3w)

# Calculate feature-specific contrast levels based on quantiles of non-zero expression.
# Particularly useful when plotting multiple markers
FeaturePlot(mouse_3w, features = c("HK1", "HK2"), 
            min.cutoff = "q10", max.cutoff = "q90", 
            cols = c("lightgrey", "red"))





