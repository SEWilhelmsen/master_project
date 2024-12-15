# Pearson's correlation coefficient (ikke ferdig)
# Input: data frame or matrix of normalized data per time point, per condition
# Output: heatmap of correlation coefficient
# Help: Jason's article https://www.sciencedirect.com/science/article/pii/S0378427423010950?via%3Dihub#fig0010
# Silje Wilhelmsen

# Load libraries
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# 1.	Perform PCA per time point and per condition.
# Load data for one time point
load(mouse_3w)

# Separate conditions (before filtering for cell types)
# Assign the column "orig.ident" as the identity columns for the Seurat object.
Idents(mouse_vcm_3w) <- mouse_vcm_3w@meta.data$orig.ident


# Do PCA
mouse_3w <- RunPCA(mouse_3w, npcs = 20)

# 2.	Save a plot of PCA (x-axis = PC1, y-axis = PC2)
# 3.	Perform Pearson correlation coefficient analysis
# 4.	Create a heatmap of the clusters per time point. 





# Assume you have already loaded your Seurat object (mouse_3w)
# Normalize the data
mouse_3w <- NormalizeData(mouse_3w)

# Identify the 2000 most variable genes (adjust as needed)
mouse_3w <- FindVariableFeatures(mouse_3w, selection.method = "vst", nfeatures = 2000)

# Scale the data
mouse_3w <- ScaleData(mouse_3w)

# Perform PCA
mouse_3w <- RunPCA(mouse_3w, features = VariableFeatures(object = mouse_3w))

# You can visualize the PCA results to decide on the number of PCs to use.
ElbowPlot(mouse_3w, ndims = 20, reduction = "pca")

# Assuming we decide to use the first 10 PCs based on the ElbowPlot
# Identify clusters in the PCA space
mouse_3w <- FindNeighbors(mouse_3w, dims = 1:5)
mouse_3w <- FindClusters(mouse_3w, resolution = 0.1)  # Adjust resolution for different clustering granularity

# Visualize the clusters (optional)
DimPlot(mouse_3w, reduction = "pca")

# Extract the PCA embeddings
pca_data <- Embeddings(mouse_3w, "pca")

# Convert PCA embeddings to a matrix for correlation analysis
pca_matrix <- as.matrix(pca_data)

# Compute Pearson's correlation matrix on the PCA components
cor_pca <- cor(pca_matrix, use = "everything", method = "pearson")

# Inspect correlation matrix
print(cor_pca)

# Compute correlation matrix if not already computed
cor_matrix <- cor(pca_matrix, use = "everything", method = "pearson")

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
        row_title = "PCA Components")