# UMAP 

# After doing add_atac_layer_to_seurat_object.R



# Prepare data
SeuratObj <- RunTFIDF(SeuratObj) # Normalize data
SeuratObj <- FindTopFeatures(SeuratObj, min.cutoff = 'q0') # Find features for dimensionality reduction
SeuratObj <- RunSVD(SeuratObj) # Perform LSI (a type of PCA for count data) before UMAP
SeuratObj <- RunUMAP(SeuratObj, reduction = "lsi", dims = 1:30)

# Create plot
DimPlot(SeuratObj, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend() 
