# Run Umap 
# Input:
# Output: 
# Silje Wilhelmsen



source("C:/Users/siljeew/snRNAseq/Analysis/UMAP.R")


# UMAP 

# Perform Pearsons correlation coefficient 

# Extract the PCA embeddings
umap_data <- Embeddings(mouse_3w, "umap")

# Convert umap embeddings to a matrix for correlation analysis
umap_matrix <- as.matrix(umap_data)

# Compute Pearson's correlation matrix on the umap components
cor_umap <- cor(umap_matrix, use = "everything", method = "pearson")

# Inspect correlation matrix
print(cor_umap)

Heatmap(cor_umap,
        name = "Pearson correlation", 
        show_column_names = TRUE,
        show_row_names = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE, 
        column_title = "Correlation", 
        row_title = "UMAP")

# Compute correlation matrix if not already computed
cor_matrix <- cor(umap_matrix, use = "everything", method = "pearson")

# Create a custom color function
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Create heatmap with ComplexHeatmap
Heatmap(cor_matrix, 
        name = "Correlation",
        col = col_fun,
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_title = "Correlation Heatmap",
        row_title = "umap Components")
