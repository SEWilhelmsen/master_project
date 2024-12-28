### RNA sequencing data SHAM vs. AB - Preprocessing
## Timepoint: 12 hours
## SHAM: 'TP232', 'TP234', 'TP265'
## AB: 'TP231', 'TP245', 'TP264'
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

counts <- Read10X_h5("D:/MultiomeSeq/TP23-2/outs/filtered_feature_bc_matrix.h5")
TP232 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "12 Hours - SHAM"
)

counts <- Read10X_h5("D:/MultiomeSeq/TP23-4/outs/filtered_feature_bc_matrix.h5")
TP234 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "12 Hours - SHAM"
)

counts <- Read10X_h5("D:/MultiomeSeq/TP26-5/outs/filtered_feature_bc_matrix.h5")
TP265 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "12 Hours - SHAM"
)

counts <- Seurat::Read10X_h5("D://MultiomeSeq/TP23-1/outs/filtered_feature_bc_matrix.h5")
TP231 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "12 Hours - AB"
)

counts <- Read10X_h5("D://MultiomeSeq/TP24-5/outs/filtered_feature_bc_matrix.h5")
TP245 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "12 Hours - AB"
)

counts <- Read10X_h5("D://MultiomeSeq/TP26-4/outs/filtered_feature_bc_matrix.h5")
TP264 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "12 Hours - AB"
)


options(future.globals.maxSize = 10 * 1024^3) # This sets the maximum size to 5 GB. 
options(timeout = 1200)

## Annotation of nuclei - Assignment to each nuclei a score on what type of cell type the nucleus comes from

Prediction_TP232 <- RunAzimuth(TP232, reference = "heartref")
Prediction_TP234 <- RunAzimuth(TP234, reference = "heartref")
Prediction_TP265 <- RunAzimuth(TP265, reference = "heartref")
Prediction_TP231 <- RunAzimuth(TP231, reference = "heartref")
Prediction_TP245 <- RunAzimuth(TP245, reference = "heartref")
Prediction_TP264 <- RunAzimuth(TP264, reference = "heartref")


## Change names on the rows of the Predictions seurat objects to nuclei ID + Animal ID#
rownames_TP232 <- rownames(Prediction_TP232@meta.data)
rownames_TP234 <- rownames(Prediction_TP234@meta.data)
rownames_TP265 <- rownames(Prediction_TP265@meta.data)
rownames_TP231 <- rownames(Prediction_TP231@meta.data)
rownames_TP245 <- rownames(Prediction_TP245@meta.data)
rownames_TP264 <- rownames(Prediction_TP264@meta.data)

new_rownames_TP232 <- gsub("-1$", "-TP232", rownames_TP232)
new_rownames_TP234 <- gsub("-1$", "-TP234", rownames_TP234)
new_rownames_TP265 <- gsub("-1$", "-TP265", rownames_TP265)
new_rownames_TP231 <- gsub("-1$", "-TP231", rownames_TP231)
new_rownames_TP245 <- gsub("-1$", "-TP245", rownames_TP245)
new_rownames_TP264 <- gsub("-1$", "-TP264", rownames_TP264)


colnames(Prediction_TP232) <- new_rownames_TP232
colnames(Prediction_TP234) <- new_rownames_TP234
colnames(Prediction_TP265) <- new_rownames_TP265
colnames(Prediction_TP231) <- new_rownames_TP231
colnames(Prediction_TP245) <- new_rownames_TP245
colnames(Prediction_TP264) <- new_rownames_TP264



## Merge all animals into one Seurat object

mouse_12h <- merge(Prediction_TP232, 
                  y=c(Prediction_TP234,
                      Prediction_TP265,
                      Prediction_TP231, 
                      Prediction_TP245,
                      Prediction_TP264), 
                  project = "mouse_12h_full_merge_project")
mouse_12h <- NormalizeData(mouse_12h, normalization.method = "LogNormalize", scale.factor = 10000)
mouse_12h <- FindVariableFeatures(mouse_12h)
mouse_12h <- ScaleData(mouse_12h)
mouse_12h <- RunPCA(mouse_12h, npcs = 20)
mouse_12h <- RunUMAP(mouse_12h, dims = 1:20)
mouse_12h <- IntegrateLayers(object = mouse_12h, 
                              method = RPCAIntegration, 
                              orig.reduction = "pca", 
                              new.reduction = "integrated.rpca",
                              verbose = FALSE)
mouse_12h <- JoinLayers(mouse_12h)


# Save Seurat object
SaveH5Seurat(mouse_12h, 
             filename = "mouse_12h.h5seurat", 
             overwrite = TRUE)
saveRDS(mouse_12h, file = "C:/Users/Labuser/snRNAseq/tmp/mouse_12h.Rds")

## Filter for ventricular cardiomyocytes
mouse_12h_vcm <- subset(mouse_12h, 
                       predicted.celltype.l2.score >= 0.7 & 
                         mapping.score >= 0.7 & 
                         predicted.celltype.l1.score >= 0.7 & 
                         mapping.score >= 0.7 &
                         predicted.celltype.l2 == 'Ventricular Cardiomycoyte')
head(mouse_12h_vcm, 5) #control post: inspect the top 5 rows in the dataset
tail(mouse_12h_vcm, 5)

SaveH5Seurat(mouse_12h_vcm, 
             filename = "mouse_12h_vcm.h5seurat", 
             overwrite = TRUE) 
saveRDS(mouse_12h_vcm, file = "C:/Users/Labuser/snRNAseq/tmp/mouse_12h_vcm.Rds")

## Find top genes
# Assign the column "orig.ident" as the identity columns for the Seurat object.
Idents(mouse_12h_vcm) <- mouse_12h_vcm@meta.data$orig.ident

# Find differentially expressed genes
mouse_12h_vcm_markers <- FindMarkers(mouse_12h_vcm, ident.1 = "12 Hours - AB", ident.2 = "12 Hours - SHAM")

mouse_12h_vcm_markers_all_genes <- mouse_12h_vcm_markers %>%
  tibble::rownames_to_column(var = "gene") %>%
  mutate(log_pval = -log10(p_val_adj))
write.xlsx(mouse_12h_vcm_markers_all_genes, file = "mouse_12h_vcm_markers_all_genes.xlsx")

# Filter for significant genes
mouse_12h_vcm_markers_significant_genes <- mouse_12h_vcm_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
write.xlsx(mouse_12h_vcm_markers_significant_genes, file = "mouse_12h_vcm_markers_significant_genes.xlsx") # Save data frame as excel
