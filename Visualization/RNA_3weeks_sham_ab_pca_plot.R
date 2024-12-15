## PCA plot
# Input:
# Output: 
# Silje Wilhelmsen

library(ggplot2)

# # To load Seurat object from .h5seurat file
# mouse_3w_sham_join_layers <- LoadH5Seurat("mouse_3w_join_layers.h5seurat")

# Plot of clusters is C8glei. Should the data be more specific? 
DimPlot(mouse_3w_pca, reduction = "pca")


# Plot the first two principal components
DimPlot(mouse_3w_pca, reduction = "pca", group.by = "orig.ident") +
  ggtitle("PCA Plot of 3 Weeks Data by Condition") +
  theme_minimal()

# Is SCT-transform different? 
mouse_3w_pca_sct <- SCTransform(mouse_3w_pca, verbose = FALSE)

DimPlot(mouse_3w_pca_sct, reduction = "pca", group.by = "orig.ident") +
  ggtitle("PCA Plot of 3 Weeks Data by Condition") +
  theme_minimal()


# DimPlot(mouse_3w_umap, reduction = "umap", group.by = "orig.ident") +
#   ggtitle("PCA Plot of 3 Weeks Data by Condition") +
#   theme_minimal()
# 
# # Save plot ???
# ggsave("Plots/dimplot.png", plot = DimPlot)