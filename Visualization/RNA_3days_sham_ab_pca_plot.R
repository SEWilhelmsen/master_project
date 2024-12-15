## PCA plot (ikke ferdig)
# Input:
# Output: 
# Silje Wilhelmsen

library(ggplot2)
library(Seurat)
library(dplyr)

# # To load Seurat object from .h5seurat file
# mouse_3d_sham_join_layers <- LoadH5Seurat("mouse_3d_join_layers.h5seurat")


# Plot the first two principal components
DimPlot(mouse_3d_pca, reduction = "pca", group.by = "orig.ident") +
  ggtitle("PCA Plot of 1 Week Data by Condition") +
  theme_minimal()

# mouse_3d_pca_sct <- SCTransform(mouse_3d_pca, verbose = FALSE)
# 
# DimPlot(mouse_3d_pca_sct, reduction = "pca", group.by = "orig.ident") +
#   ggtitle("PCA Plot of 1 Week Data by Condition") +
#   theme_minimal()
# 
# 
# DimPlot(mouse_3d_umap, reduction = "umap", group.by = "orig.ident") +
#   ggtitle("PCA Plot of 1 Week1 Data by Condition") +
#   theme_minimal()


# Extract the standard deviations of the principal components
stdevs <- mouse_3d_pca[["pca"]]@stdev

# Calculate the percentage of variance explained by each PC
explained_variance <- (stdevs^2) / sum(stdevs^2) * 100

# Create a data frame for easier handling and plotting
explained_variance_df <- data.frame(
  PC = paste0("PC", 1:length(explained_variance)),
  Variance = explained_variance
)

# Sort the data frame in descending order of explained variance
explained_variance_df <- explained_variance_df[order(-explained_variance_df$Variance),]


# Inspect the first few rows of the explained variance data frame
head(explained_variance_df)



# Define custom colors (modify as needed based on number of PCs)
custom_colors <- c("PC1" = "#1f77b4", "PC2" = "#ff7f0e", "PC3" = "#2ca02c", "PC4" = "#d62728",
                   "PC5" = "#9467bd", "PC6" = "#8c564b", "PC7" = "#e377c2", "PC8" = "#7f7f7f",
                   "PC9" = "#bcbd22", "PC10" = "#17becf", "PC11" = "#aec7e8", "PC12" = "#ffbb78",
                   "PC13" = "#98df8a", "PC14" = "#ff9896", "PC15" = "#c5b0d5", "PC16" = "#c49c94",
                   "PC17" = "#f7b6d2", "PC18" = "#c7c7c7", "PC19" = "#dbdb8d", "PC20" = "#9edae5")

# Plot the explained variance with custom colors
variance_plot <- ggplot(explained_variance_df, aes(x = reorder(PC, -Variance), y = Variance, fill = PC)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(title = "Explained Variance by Principal Components 3 Days",
       x = "Principal Components",
       y = "Percentage of Variance Explained") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Display the plot
print(variance_plot)

# # Extract PCA loadings
# mouse_3d_pca_loadings <- Loadings(mouse_3d_pca, reduction = "pca")
# mouse_3d_pca_loadings_df <- as.data.frame(mouse_3d_pca_loadings)
# 
# # Save complete PCA loadings to CSV file
# write.csv(mouse_3d_pca_loadings_df, "mouse_3d_pca_loadings_df.csv", 
#           row.names = TRUE)
