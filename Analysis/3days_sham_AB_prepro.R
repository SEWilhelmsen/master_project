### RNA sequencing data SHAM vs. AB - Preprocessing 
## Timepoint: 3 days
## SHAM: 'TP9-1', 'TP9-5', 'TP11-5'
## AB: 'TP9-2', 'TP10-2', 'TP11-1A'
## Silje Wilhelmsen


## Relevant litterature https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.117.030742

library(SoupX)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(biovizBase)
library(Seurat)
library(Azimuth)
library(SeuratData)
library(TFBSTools)
library(spacexr)
library(Matrix)
library(zellkonverter)
library(SingleCellExperiment)
library(reticulate)
library(harmony)
library(SeuratDisk)
library(somepackage)
library(clusterProfiler)
library(org.Hs.eg.db)
library(heartref.SeuratData)
library(writexl)

options(future.globals.maxSize = 9 * 1024^3)  # Setting max size to 2 GiB in case of error due to size.

## Load the data - 3 sham and 3 AB mice - 3 days after ORAB operation

counts <- Seurat::Read10X_h5("W:/CardioTarget/Multiome seq data/TP9-1/outs/filtered_feature_bc_matrix.h5")
TP91 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "3 Days - SHAM"
)

counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/TP9-5/outs/filtered_feature_bc_matrix.h5")
TP95 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "3 Days - SHAM"
)

counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/TP11-5/outs/filtered_feature_bc_matrix.h5")
TP115 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "3 Days - SHAM"
)

counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/TP9-2/outs/filtered_feature_bc_matrix.h5")
TP92 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "3 Days - AB"
)

counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/TP10-2/outs/filtered_feature_bc_matrix.h5")
TP102 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "3 Days - AB"
)

counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/TP11-1A/outs/filtered_feature_bc_matrix.h5")
TP111 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "3 Days - AB"
)


## Annotation of nuclei - Assignment to each nuclei a score on what type of cell type the nucleus comes from
Prediction_TP91 <- RunAzimuth(TP91, reference = "heartref")
Prediction_TP95 <- RunAzimuth(TP95, reference = "heartref")
Prediction_TP115 <- RunAzimuth(TP115, reference = "heartref")
Prediction_TP92 <- RunAzimuth(TP92, reference = "heartref")
Prediction_TP102 <- RunAzimuth(TP102, reference = "heartref")
Prediction_TP111 <- RunAzimuth(TP111, reference = "heartref")


## Change names on the rows of the Predictions seurat objects to nuclei ID + Animal ID#
rownames_TP91 <- rownames(Prediction_TP91@meta.data)
rownames_TP95 <- rownames(Prediction_TP95@meta.data)
rownames_TP115 <- rownames(Prediction_TP115@meta.data)
rownames_TP92 <- rownames(Prediction_TP92@meta.data)
rownames_TP102 <- rownames(Prediction_TP102@meta.data)
rownames_TP111 <- rownames(Prediction_TP111@meta.data)

new_rownames_TP91 <- gsub("-1$", "-TP91", rownames_TP91)
new_rownames_TP95 <- gsub("-1$", "-TP95", rownames_TP95)
new_rownames_TP115 <- gsub("-1$", "-TP115", rownames_TP115)
new_rownames_TP92 <- gsub("-1$", "-TP92", rownames_TP92)
new_rownames_TP102 <- gsub("-1$", "-TP102", rownames_TP102)
new_rownames_TP111 <- gsub("-1$", "-TP111", rownames_TP111) 

colnames(Prediction_TP91) <- new_rownames_TP91
colnames(Prediction_TP95) <- new_rownames_TP95
colnames(Prediction_TP115) <- new_rownames_TP115
colnames(Prediction_TP92) <- new_rownames_TP92
colnames(Prediction_TP102) <- new_rownames_TP102
colnames(Prediction_TP111) <- new_rownames_TP111

## Merge all 6 animals into one Seurat object

mouse_3d <- merge(Prediction_TP91, 
                  y=c(Prediction_TP95, 
                      Prediction_TP115, 
                      Prediction_TP92, 
                      Prediction_TP102, 
                      Prediction_TP111), 
                  project = "mouse_3d_Full_Merge_Project")
mouse_3d_normalize <- NormalizeData(mouse_3d, 
                                    normalization.method = "LogNormalize", 
                                    scale.factor = 10000)
mouse_3d_find_features <- FindVariableFeatures(mouse_3d_normalize)
mouse_3d_scale <- ScaleData(mouse_3d_find_features)
mouse_3d_pca <- RunPCA(mouse_3d_scale, npcs = 20)
mouse_3d_umap <- RunUMAP(mouse_3d_pca, dims = 1:20)

mouse_3d_integrate <- IntegrateLayers(object = mouse_3d_umap, 
                            method = RPCAIntegration, 
                            orig.reduction = "pca", 
                            new.reduction = "integrated.rpca",
                            verbose = FALSE)
mouse_3d_join_layers <- JoinLayers(mouse_3d_integrate)


# Save Seurat object
SaveH5Seurat(mouse_3d_join_layers, 
             file = "mouse_3d_join_layers.h5seurat", 
             overwrite = TRUE, 
             verbose = TRUE)


## Filter for ventricular cardiomyocytes?
mouse_3d_vcm <- subset(mouse_3d_join_layers, 
                          predicted.celltype.l2.score >= 0.7 & 
                            mapping.score >= 0.7 & 
                            predicted.celltype.l1.score >= 0.7 & 
                            mapping.score >= 0.7 &
                            predicted.celltype.l2 == 'Ventricular Cardiomycoyte')
head(mouse_3d_vcm, 5) #control post: inspect the top 5 rows in the dataset
tail(mouse_3d_vcm, 5)

# Save Seurat object
SaveH5Seurat(mouse_3d_vcm, 
             filename = "mouse_3d_vcm.h5seurat", 
             overwrite = TRUE) 


# Assign the column "orig.ident" as the identity columns for the Seurat object.
Idents(mouse_3d_vcm) <- mouse_3d_vcm@meta.data$orig.ident

# Find differentially expressed genes
mouse_3d_vcm_markers <- FindMarkers(mouse_3d_vcm, ident.1 = "3 Days - AB", ident.2 = "3 Days - SHAM")

mouse_3d_vcm_markers_all_genes <- mouse_3d_vcm_markers %>%
  tibble::rownames_to_column(var = "gene") %>%
  mutate(log_pval = -log10(p_val_adj))
write.xlsx(mouse_3d_vcm_markers_all_genes, file = "mouse_3d_vcm_markers_all_genes.xlsx")

# Filter for significant genes
mouse_3d_vcm_markers_significant_genes <- mouse_3d_vcm_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
write.xlsx(mouse_3d_vcm_markers_significant_genes, file = "mouse_3d_vcm_markers_significant_genes.xlsx") # Save data frame as excel

# Filter for top 200 significant genes
mouse_3d_vcm_markers_significant_genes_top_200 <- mouse_3d_vcm_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(200)
write.xlsx(mouse_3d_vcm_markers_significant_genes_top_200, file = "mouse_3d_vcm_markers_significant_genes_top_200.xlsx") # Save data frame as excel

