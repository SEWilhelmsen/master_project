### RNA sequencing data SHAM vs. AB - Preprocessing
## Timepoint: 3 weeks
## SHAM: 'TP5-1', 'Sample1-TP75B', 'Sample3-TP52A'
## AB: 'Sample2-TP511B', Sample4-TP62A', 'TP51-3'
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
library(openxlsx)



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

counts <- Seurat::Read10X_h5("W:/CardioTarget/Multiome seq data/GCF0742_GCF0795_CellRanger/Sample2-TP511B/outs/filtered_feature_bc_matrix.h5")
TP511 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "3 Weeks - AB"
)

counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/TP51-3/outs/filtered_feature_bc_matrix.h5")
TP513 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "3 Weeks - AB"
)


counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/GCF0742_GCF0795_CellRanger/Sample4-TP62A/outs/filtered_feature_bc_matrix.h5")
TP62 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "3 Weeks - AB"
)

options(future.globals.maxSize = 5 * 1024^3) # This sets the maximum size to 5 GB. 

## Annotation of nuclei - Assignment to each nuclei a score on what type of cell type the nucleus comes from

Prediction_TP51 <- RunAzimuth(TP51, reference = "heartref")
Prediction_TP75 <- RunAzimuth(TP75, reference = "heartref")
Prediction_TP52 <- RunAzimuth(TP52, reference = "heartref")
Prediction_TP511 <- RunAzimuth(TP511, reference = "heartref")
Prediction_TP62 <- RunAzimuth(TP62, reference = "heartref")
Prediction_TP513 <- RunAzimuth(TP513, reference = "heartref")


## Change names on the rows of the Predictions seurat objects to nuclei ID + Animal ID#
rownames_TP51 <- rownames(Prediction_TP51@meta.data)
rownames_TP75 <- rownames(Prediction_TP75@meta.data)
rownames_TP52 <- rownames(Prediction_TP52@meta.data)
rownames_TP511 <- rownames(Prediction_TP511@meta.data)
rownames_TP62 <- rownames(Prediction_TP62@meta.data)
rownames_TP513 <- rownames(Prediction_TP513@meta.data)

new_rownames_TP51 <- gsub("-1$", "-TP51", rownames_TP51)
new_rownames_TP75 <- gsub("-1$", "-TP75", rownames_TP75)
new_rownames_TP52 <- gsub("-1$", "-TP52", rownames_TP52)
new_rownames_TP511 <- gsub("-1$", "-TP511", rownames_TP511)
new_rownames_TP62 <- gsub("-1$", "-TP62", rownames_TP62)
new_rownames_TP513 <- gsub("-1$", "-TP513", rownames_TP513)

colnames(Prediction_TP51) <- new_rownames_TP51
colnames(Prediction_TP75) <- new_rownames_TP75
colnames(Prediction_TP52) <- new_rownames_TP52
colnames(Prediction_TP511) <- new_rownames_TP511
colnames(Prediction_TP62) <- new_rownames_TP62
colnames(Prediction_TP513) <- new_rownames_TP513

## Merge all 6 animals into one Seurat object

mouse_3w <- merge(Prediction_TP51, 
                  y=c(Prediction_TP75, 
                      Prediction_TP52, 
                      Prediction_TP511, 
                      Prediction_TP62, 
                      Prediction_TP513), 
                  project = "mouse_3w_full_merge_project")
mouse_3w <- NormalizeData(mouse_3w, 
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000)
mouse_3w <- FindVariableFeatures(mouse_3w)
mouse_3w <- ScaleData(mouse_3w)
mouse_3w <- RunPCA(mouse_3w, npcs = 20)
mouse_3w <- RunUMAP(mouse_3w, dims = 1:20)
mouse_3w <- IntegrateLayers(object = mouse_3w, 
                            method = RPCAIntegration, 
                            orig.reduction = "pca", 
                            new.reduction = "integrated.rpca",
                            verbose = FALSE)
mouse_3w_join_layers <- JoinLayers(mouse_3w)


# Save Seurat object
SaveH5Seurat(mouse_3w_join_layers, 
             filename = "mouse_3w_join_layers.h5seurat", 
             overwrite = TRUE)


## Filter for ventricular cardiomyocytes
mouse_3w_vcm <- subset(mouse_3w_join_layers, 
                       predicted.celltype.l2.score >= 0.7 & 
                         mapping.score >= 0.7 & 
                         predicted.celltype.l1.score >= 0.7 & 
                         mapping.score >= 0.7 &
                         predicted.celltype.l2 == 'Ventricular Cardiomycoyte')
head(mouse_3w_vcm, 5) #control post: inspect the top 5 rows in the dataset
tail(mouse_3w_vcm, 5)

SaveH5Seurat(mouse_3w_vcm, 
             filename = "mouse_3w_vcm.h5seurat", 
             overwrite = TRUE) 

## Find top genes
# Assign the column "orig.ident" as the identity columns for the Seurat object.
Idents(mouse_3w_vcm) <- mouse_3w_vcm@meta.data$orig.ident

# Find differentially expressed genes
mouse_3w_vcm_markers <- FindMarkers(mouse_3w_vcm, ident.1 = "3 Weeks - AB", ident.2 = "3 Weeks - SHAM")

mouse_3w_vcm_markers_all_genes <- mouse_3w_vcm_markers %>%
  tibble::rownames_to_column(var = "gene") %>%
  mutate(log_pval = -log10(p_val_adj))
write.xlsx(mouse_3w_vcm_markers_all_genes, file = "mouse_3w_vcm_markers_all_genes.xlsx")

# Filter for significant genes
mouse_3w_vcm_markers_significant_genes <- mouse_3w_vcm_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
write.xlsx(mouse_3w_vcm_markers_significant_genes, file = "mouse_3w_vcm_markers_significant_genes.xlsx") # Save data frame as excel

# # Filter for top 200 significant genes
# mouse_3w_vcm_markers_significant_genes_top_200 <- mouse_3w_vcm_markers %>%
#   dplyr::filter(p_val_adj < 0.05) %>%
#   arrange(desc(abs(avg_log2FC))) %>%
#   head(200)
# write.xlsx(mouse_3w_vcm_markers_significant_genes_top_200, file = "mouse_3w_vcm_markers_significant_genes_top_200.xlsx") # Save data frame as excel
