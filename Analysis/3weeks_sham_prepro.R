### RNA sequencing data SHAM - Preprocessing (ikke ferdig)
## Timepoint: 3 weeks
## SHAM: 'TP5-1', 'Sample1-TP75B', 'Sample3-TP52A'
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



## Load the data - 3 sham and 3 AB mice - 3 weeks after ORAB operation

counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/TP5-1/outs/filtered_feature_bc_matrix.h5")
TP51 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "3 Weeks - SHAM"
)

counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/GCF0742_GCF0795_CellRanger/Sample1-TP75B/outs/filtered_feature_bc_matrix.h5")
TP75 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "3 Weeks - SHAM"
)

counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/GCF0742_GCF0795_CellRanger/Sample3-TP52A/outs/filtered_feature_bc_matrix.h5")
TP52 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "3 Weeks - SHAM"
)


## Annotation of nuclei - Assignment to each nuclei a score on what type of cell type the nucleus comes from

Prediction_TP51 <- RunAzimuth(TP51, reference = "heartref")
Prediction_TP75 <- RunAzimuth(TP75, reference = "heartref")
Prediction_TP52 <- RunAzimuth(TP52, reference = "heartref")

## Change names on the rows of the Predictions seurat objects to nuclei ID + Animal ID#
rownames_TP51 <- rownames(Prediction_TP51@meta.data)
rownames_TP75 <- rownames(Prediction_TP75@meta.data)
rownames_TP52 <- rownames(Prediction_TP52@meta.data)

new_rownames_TP51 <- gsub("-1$", "-TP51", rownames_TP51)
new_rownames_TP75 <- gsub("-1$", "-TP75", rownames_TP75)
new_rownames_TP52 <- gsub("-1$", "-TP52", rownames_TP52)

colnames(Prediction_TP51) <- new_rownames_TP51
colnames(Prediction_TP75) <- new_rownames_TP75
colnames(Prediction_TP52) <- new_rownames_TP52

## Merge all 6 animals into one Seurat object

mouse_3w_sham <- merge(Prediction_TP51, 
                  y=c(Prediction_TP75, 
                      Prediction_TP52), 
                  project = "Full_Merge_Project_sham_3w")
mouse_3w_sham_normalize <- NormalizeData(mouse_3w_sham, 
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000)
## If you want to save normalized data
# Extract normalized data matrix
mouse_3w_sham_normalize <- GetAssayData(mouse_3w_sham_normalize, assay = "RNA", slot = "data")
mouse_3w_sham_normalize_df <- as.data.frame(as.matrix(mouse_3w_sham_normalize)) # Convert to data frame for saving as CSV
write.csv(mouse_3w_sham_normalize_df, "mouse_3w_sham_normalize_df.csv", row.names = TRUE) # Save as CSV

# Continue

mouse_3w_sham_find_features <- FindVariableFeatures(mouse_3w_sham_normalize)
mouse_3w_sham_scale <- ScaleData(mouse_3w_sham_find_features)
mouse_3w_sham_pca <- RunPCA(mouse_3w_sham_scale, npcs = 20)
mouse_3w_sham_umap <- RunUMAP(mouse_3w_sham_pca, dims = 1:20)
mouse_3w_sham_integrate <- IntegrateLayers(object = mouse_3w_sham_umap, 
                            method = RPCAIntegration, 
                            orig.reduction = "pca", 
                            new.reduction = "integrated.rpca",
                            verbose = FALSE)
mouse_3w_sham_join_layers <- JoinLayers(mouse_3w_sham_integrate)


# Save Seurat object
SaveH5Seurat(mouse_3w_sham_join_layers, overwrite = TRUE)

mouse_3w_sham_join_layers <- GetAssayData(mouse_3w_sham_join_layers, assay = "RNA", slot = "data") # "data" slot contains the normalized data
mouse_3w_sham_matrix <- as.matrix(mouse_3w_sham_join_layers)

## Filter for ventricular cardiomyocytes
mouse_3w_sham_vcm <- subset(mouse_3w_sham, 
                       predicted.celltype.l2.score >= 0.7 & 
                         mapping.score >= 0.7 & 
                         predicted.celltype.l1.score >= 0.7 & 
                         mapping.score >= 0.7 &
                         predicted.celltype.l2 == 'Ventricular Cardiomycoyte')
head(mouse_3w_sham_vcm, 5) #control post: inspect the top 5 rows in the dataset
tail(mouse_3w_sham_vcm, 5)

SaveH5Seurat(mouse_3w_sham_vcm, "V:/Silje/snRNAseq/3weeks_sham_ab/mouse_3w_sham_vcm", overwrite = TRUE) # Save Seurat object?




## Find top genes
# Assign the column "orig.ident" as the identity columns for the Seurat object.
Idents(mouse_3w_sham_vcm) <- mouse_3w_sham_vcm@meta.data$orig.ident

# Find differentially expressed genes
mouse_3w_sham_vcm_markers <- FindMarkers(mouse_3w_sham_vcm, ident.1 = "3 Weeks - AB", ident.2 = "3 Weeks - SHAM")

mouse_3w_sham_vcm_markers <- markers %>%
  tibble::rownames_to_column(var = "gene") %>%
  mutate(log_pval = -log10(p_val_adj))

mouse_3w_sham_vcm_markers_significant_genes <- markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))

write_xlsx(mouse_3w_sham_vcm_markers_significant_genes, "V:/Silje/snRNAseq/3weeks_sham_ab/mouse_3w_sham_vcm_markers_significant_genes_3w.xlsx") # Save data frame as excel

# mouse_3w_sham_vcm_markers_top_genes <- markers %>%
#   dplyr::filter(p_val_adj < 0.05) %>%
#   arrange(desc(abs(avg_log2FC))) %>%
#   head(200)

