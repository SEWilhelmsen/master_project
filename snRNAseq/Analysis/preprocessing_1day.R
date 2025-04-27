### RNA sequencing data SHAM vs. ORAB - Preprocessing
## Timepoint: 1 Day
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
library(openxlsx)


## Load the data after ORAB operation

counts <- Seurat::Read10X_h5("D:/MultiomeSeq/TP17-2/outs/filtered_feature_bc_matrix.h5")
TP172 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "1 Day - SHAM"
)


counts <- Seurat::Read10X_h5("D:/MultiomeSeq/TP13-1/outs/filtered_feature_bc_matrix.h5")
TP131 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "1 Day - SHAM"
)

counts <- Seurat::Read10X_h5("D:/MultiomeSeq/TP17-3/outs/filtered_feature_bc_matrix.h5")
TP173 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "1 Day - SHAM"
)


counts <- Read10X_h5("D:/MultiomeSeq/TP14-2/outs/filtered_feature_bc_matrix.h5")
TP142 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "1 Day - ORAB"
)

counts <- Read10X_h5("D:/MultiomeSeq/TP14-1/outs/filtered_feature_bc_matrix.h5")
TP141 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "1 Day - ORAB"
)

counts <- Read10X_h5("D:/MultiomeSeq/TP15-3/outs/filtered_feature_bc_matrix.h5")
TP153 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "1 Day - ORAB"
)


## Annotation of nuclei - Assignment to each nuclei a score on what type of cell type the nucleus comes from

options(future.globals.maxSize = 10 * 1024^3) # If error on the next lines, change max size to 10 GB. 
getOption("future.globals.maxSize")

Prediction_TP172 <- RunAzimuth(TP172, reference = "heartref")
Prediction_TP131 <- RunAzimuth(TP131, reference = "heartref")
Prediction_TP173 <- RunAzimuth(TP173, reference = "heartref")
Prediction_TP142 <- RunAzimuth(TP142, reference = "heartref")
Prediction_TP141 <- RunAzimuth(TP141, reference = "heartref")
Prediction_TP153 <- RunAzimuth(TP153, reference = "heartref")


## Change names on the rows of the Predictions seurat objects to nuclei ID + Animal ID#
rownames_TP172 <- rownames(Prediction_TP172@meta.data)
rownames_TP131 <- rownames(Prediction_TP131@meta.data)
rownames_TP173 <- rownames(Prediction_TP173@meta.data)
rownames_TP142 <- rownames(Prediction_TP142@meta.data)
rownames_TP141 <- rownames(Prediction_TP141@meta.data)
rownames_TP153 <- rownames(Prediction_TP153@meta.data)

new_rownames_TP172 <- gsub("-1$", "-TP172", rownames_TP172)
new_rownames_TP131 <- gsub("-1$", "-TP131", rownames_TP131)
new_rownames_TP173 <- gsub("-1$", "-TP173", rownames_TP173)
new_rownames_TP142 <- gsub("-1$", "-TP142", rownames_TP142)
new_rownames_TP141 <- gsub("-1$", "-TP141", rownames_TP141)
new_rownames_TP153 <- gsub("-1$", "-TP153", rownames_TP153)

colnames(Prediction_TP172) <- new_rownames_TP172
colnames(Prediction_TP131) <- new_rownames_TP131
colnames(Prediction_TP173) <- new_rownames_TP173
colnames(Prediction_TP142) <- new_rownames_TP142
colnames(Prediction_TP141) <- new_rownames_TP141
colnames(Prediction_TP153) <- new_rownames_TP153

## Merge all animals into one Seurat object

mouse_1d <- merge(Prediction_TP172, 
                  y=c(Prediction_TP131,
                      Prediction_TP173,
                      Prediction_TP142, 
                      Prediction_TP141,
                      Prediction_TP153), 
                  project = "mouse_1d_full_merge_project")
mouse_1d <- NormalizeData(mouse_1d, 
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000)
mouse_1d <- FindVariableFeatures(mouse_1d)
mouse_1d <- ScaleData(mouse_1d)
mouse_1d <- RunPCA(mouse_1d, npcs = 20)
mouse_1d <- RunUMAP(mouse_1d, dims = 1:20)
mouse_1d <- IntegrateLayers(object = mouse_1d, 
                            method = RPCAIntegration, 
                            orig.reduction = "pca", 
                            new.reduction = "integrated.rpca",
                            verbose = FALSE)
mouse_1d <- JoinLayers(mouse_1d)


# Save Seurat object
SaveH5Seurat(mouse_1d, 
             filename = "mouse_1d.h5seurat", 
             overwrite = TRUE)
saveRDS(mouse_1d, file = "C:/Users/Labuser/snRNAseq/tmp/mouse_1d.Rds")

## Filter for ventricular cardiomyocytes
mouse_1d_vcm <- subset(mouse_1d, 
                       predicted.celltype.l2.score >= 0.7 & 
                         mapping.score >= 0.7 & 
                         predicted.celltype.l1.score >= 0.7 & 
                         mapping.score >= 0.7 &
                         predicted.celltype.l2 == 'Ventricular Cardiomycoyte')
head(mouse_1d_vcm, 5) #control post: inspect the top 5 rows in the dataset
tail(mouse_1d_vcm, 5)

SaveH5Seurat(mouse_1d_vcm, 
             filename = "mouse_1d_vcm.h5seurat", 
             overwrite = TRUE) 
saveRDS(mouse_1d_vcm, file = "C:/Users/Labuser/snRNAseq/tmp/mouse_1d_vcm.Rds")

## Find top genes
# Assign the column "orig.ident" as the identity columns for the Seurat object.
Idents(mouse_1d_vcm) <- mouse_1d_vcm@meta.data$orig.ident

# Find differentially expressed genes
mouse_1d_vcm_markers <- FindMarkers(mouse_1d_vcm, ident.1 = "1 Day - ORAB", ident.2 = "1 Day - SHAM", logfc.threshold = 0)

mouse_1d_vcm_markers_all_genes <- mouse_1d_vcm_markers %>%
  tibble::rownames_to_column(var = "gene") %>%
  mutate(log_pval = -log10(p_val_adj))
write.xlsx(mouse_1d_vcm_markers_all_genes, file = "mouse_1d_vcm_markers_all_genes.xlsx")


# Filter for significant genes
mouse_1d_vcm_markers_significant_genes <- mouse_1d_vcm_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
write.xlsx(mouse_1d_vcm_markers_significant_genes, file = "mouse_1d_vcm_markers_significant_genes.xlsx")