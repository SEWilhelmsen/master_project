## PCA plot
# Input:
# Output: 
# Silje Wilhelmsen

library(ggplot2)

# # To load Seurat object from .h5seurat file
# mouse_1w_sham_join_layers <- LoadH5Seurat("mouse_1w_join_layers.h5seurat")


# Plot the first two principal components
DimPlot(mouse_1w_pca, reduction = "pca", group.by = "orig.ident") +
  ggtitle("PCA Plot of 1 Week Data by Condition") +
  theme_minimal()

mouse_1w_pca_sct <- SCTransform(mouse_3w_pca, verbose = FALSE)

DimPlot(mouse_1w_pca_sct, reduction = "pca", group.by = "orig.ident") +
  ggtitle("PCA Plot of 1 Week Data by Condition") +
  theme_minimal()


DimPlot(mouse_1w_umap, reduction = "umap", group.by = "orig.ident") +
  ggtitle("PCA Plot of 1 Week1 Data by Condition") +
  theme_minimal()

