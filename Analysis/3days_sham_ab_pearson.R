## Pearson correlation coefficient (ikke ferdig)
# Input:
# Output: 
# Silje Wilhelmsen


# Load necessary libraries
library(Seurat)
library(ComplexHeatmap)
library(circlize)


# # Here you should save the file? or do a Pearson correlation analysis
# mouse_3d_pca_corr <- cor(mouse_3d_pca, method = "pearson")


# Assuming mouse_3d_pca is your Seurat object after running PCA

# Extract PCA embeddings
pca_embeddings_umap <- Embeddings(mouse_3d_umap, reduction = "umap")

# Calculate the Pearson correlation matrix for the PCA embeddings
pca_cor_matrix <- cor(pca_embeddings_umap, method = "pearson")

# Print the first few rows and columns of the correlation matrix for inspection
print(head(pca_cor_matrix))

# Create a custom color function for the heatmap
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Create and display the heatmap using ComplexHeatmap
Heatmap(pca_cor_matrix,
        name = "Correlation",
        col = col_fun,
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_title = "Pearson Correlation of PCA Components",
        row_title = "PCA Components"
)

##################### scatterplot

# Load necessary libraries
library(Seurat)
library(ggplot2)

# Assuming mouse_3d_pca is your Seurat object after running PCA

# Extract PCA embeddings
pca_embeddings <- Embeddings(mouse_3d_pca, reduction = "pca")
pca_df <- as.data.frame(pca_embeddings)

# Add metadata to the data frame (e.g., conditions)
pca_df$condition <- mouse_3d_pca@meta.data$orig.ident

# Create scatter plot for the first two principal components using ggplot2
pca_scatter_plot <- ggplot(pca_df, aes(x = PC_1, y = PC_2, color = condition)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA Plot of Mouse Data", x = "PC1", y = "PC2", color = "Condition") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the scatter plot
print(pca_scatter_plot)