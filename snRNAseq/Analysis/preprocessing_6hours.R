### RNA sequencing data SHAM vs. AB - Preprocessing
## Timepoint: 6 hours
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



## Load the data 

counts <- Read10X_h5("D:/MultiomeSeq/TP22-2/outs/filtered_feature_bc_matrix.h5")
TP222 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "6 Hours - SHAM"
)

counts <- Read10X_h5("D:/MultiomeSeq/TP21-3/outs/filtered_feature_bc_matrix.h5")
TP213 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "6 Hours - SHAM"
)

counts <- Read10X_h5("D:/MultiomeSeq/TP22-1/outs/filtered_feature_bc_matrix.h5")
TP221 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "6 Hours - SHAM"
)

counts <- Seurat::Read10X_h5("D:/MultiomeSeq/TP19-2/outs/filtered_feature_bc_matrix.h5")
TP192 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "6 Hours - AB"
)

counts <- Read10X_h5("D:/MultiomeSeq/TP19-4/outs/filtered_feature_bc_matrix.h5")
TP194 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "6 Hours - AB"
)


counts <- Read10X_h5("D:/MultiomeSeq/TP22-3/outs/filtered_feature_bc_matrix.h5")
TP223 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "6 Hours - AB"
)

options(future.globals.maxSize = 10 * 1024^3) # This sets the maximum size to 5 GB. 
options(timeout = 1200) # Set time limit to 20 minutes

## Annotation of nuclei - Assignment to each nuclei a score on what type of cell type the nucleus comes from

Prediction_TP222 <- RunAzimuth(TP222, reference = "heartref")
Prediction_TP213 <- RunAzimuth(TP213, reference = "heartref")
Prediction_TP221 <- RunAzimuth(TP221, reference = "heartref")
Prediction_TP192 <- RunAzimuth(TP192, reference = "heartref")
Prediction_TP194 <- RunAzimuth(TP194, reference = "heartref")
Prediction_TP223 <- RunAzimuth(TP223, reference = "heartref")


## Change names on the rows of the Predictions seurat objects to nuclei ID + Animal ID#
rownames_TP222 <- rownames(Prediction_TP222@meta.data)
rownames_TP213 <- rownames(Prediction_TP213@meta.data)
rownames_TP221 <- rownames(Prediction_TP221@meta.data)
rownames_TP192 <- rownames(Prediction_TP192@meta.data)
rownames_TP194 <- rownames(Prediction_TP194@meta.data)
rownames_TP223 <- rownames(Prediction_TP223@meta.data)

new_rownames_TP222 <- gsub("-1$", "-TP222", rownames_TP222)
new_rownames_TP213 <- gsub("-1$", "-TP213", rownames_TP213)
new_rownames_TP221 <- gsub("-1$", "-TP221", rownames_TP221)
new_rownames_TP192 <- gsub("-1$", "-TP192", rownames_TP192)
new_rownames_TP194 <- gsub("-1$", "-TP194", rownames_TP194)
new_rownames_TP223 <- gsub("-1$", "-TP223", rownames_TP223)

colnames(Prediction_TP222) <- new_rownames_TP222
colnames(Prediction_TP213) <- new_rownames_TP213
colnames(Prediction_TP221) <- new_rownames_TP221
colnames(Prediction_TP192) <- new_rownames_TP192
colnames(Prediction_TP194) <- new_rownames_TP194
colnames(Prediction_TP223) <- new_rownames_TP223

## Merge all 6 animals into one Seurat object

mouse_6h <- merge(Prediction_TP222, 
                  y=c(Prediction_TP213, 
                      Prediction_TP221, 
                      Prediction_TP192, 
                      Prediction_TP194, 
                      Prediction_TP223), 
                  project = "mouse_6h_full_merge_project")
mouse_6h <- NormalizeData(mouse_6h, normalization.method = "LogNormalize", scale.factor = 10000)
mouse_6h <- FindVariableFeatures(mouse_6h)
mouse_6h <- ScaleData(mouse_6h)
mouse_6h <- RunPCA(mouse_6h, npcs = 20)
mouse_6h <- RunUMAP(mouse_6h, dims = 1:20)
mouse_6h <- IntegrateLayers(object = mouse_6h, 
                            method = RPCAIntegration, 
                            orig.reduction = "pca", 
                            new.reduction = "integrated.rpca",
                            verbose = FALSE)
mouse_6h <- JoinLayers(mouse_6h)


# Save Seurat object
SaveH5Seurat(mouse_6h, 
             filename = "mouse_6h.h5seurat", 
             overwrite = TRUE)

# Save RDS
saveRDS(mouse_6h, file = "C:/Users/Labuser/snRNAseq/tmp/mouse_6h.Rds")



## Filter for ventricular cardiomyocytes
mouse_6h_vcm <- subset(mouse_6h, 
                       predicted.celltype.l2.score >= 0.7 & 
                         mapping.score >= 0.7 & 
                         predicted.celltype.l1.score >= 0.7 & 
                         mapping.score >= 0.7 &
                         predicted.celltype.l2 == 'Ventricular Cardiomycoyte')
head(mouse_6h_vcm, 5) #control post: inspect the top 5 rows in the dataset
tail(mouse_6h_vcm, 5)

SaveH5Seurat(mouse_6h_vcm, filename = "mouse_6h_vcm.h5seurat", overwrite = TRUE) 

# Save RDS
saveRDS(mouse_6h_vcm, file = "C:/Users/Labuser/snRNAseq/tmp/mouse_6h_vcm.Rds")

## Find top genes
# Assign the column "orig.ident" as the identity columns for the Seurat object.
Idents(mouse_6h_vcm) <- mouse_6h_vcm@meta.data$orig.ident

# Find differentially expressed genes
mouse_6h_vcm_markers <- FindMarkers(mouse_6h_vcm, ident.1 = "6 Hours - AB", ident.2 = "6 Hours - SHAM", logfc.threshold = 0)

mouse_6h_vcm_markers_all_genes <- mouse_6h_vcm_markers %>%
  tibble::rownames_to_column(var = "gene") %>%
  mutate(log_pval = -log10(p_val_adj))
write.xlsx(mouse_6h_vcm_markers_all_genes, file = "mouse_6h_vcm_markers_all_genes.xlsx")

# Filter for significant genes
mouse_6h_vcm_markers_significant_genes <- mouse_6h_vcm_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
write.xlsx(mouse_6h_vcm_markers_significant_genes, file = "mouse_6h_vcm_markers_significant_genes.xlsx") # Save data frame as excel

