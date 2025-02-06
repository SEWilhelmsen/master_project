### RNA sequencing data SHAM vs. AB - Preprocessing
## Timepoint: 1 Week
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

# Set parameters
options(future.globals.maxSize = 5 * 1024^3) # This sets the maximum size to 5 GB. 

## Load the data - 3 sham and 3 AB mice - 1 week after ORAB operation

counts <- Seurat::Read10X_h5("W:/CardioTarget/Multiome seq data/TP3-1/outs/filtered_feature_bc_matrix.h5")
TP31 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "1 Week - SHAM"
)

counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/TP3-4/outs/filtered_feature_bc_matrix.h5")
TP34 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "1 Week - SHAM"
)

counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/TP4-1/outs/filtered_feature_bc_matrix.h5")
TP41 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "1 Week - SHAM"
)

counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/TP2-4/outs/filtered_feature_bc_matrix.h5")
TP24 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "1 Week - AB"
)

counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/TP52-2/outs/filtered_feature_bc_matrix.h5")
TP522 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "1 Week - AB"
)

counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/TP58-4/outs/filtered_feature_bc_matrix.h5")
TP584 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "1 Week - AB"
)



## Annotation of nuclei - Assignment to each nuclei a score on what type of cell type the nucleus comes from

Prediction_TP31 <- RunAzimuth(TP31, reference = "heartref")
Prediction_TP34 <- RunAzimuth(TP34, reference = "heartref")
Prediction_TP41 <- RunAzimuth(TP41, reference = "heartref")
Prediction_TP24 <- RunAzimuth(TP24, reference = "heartref")
Prediction_TP522 <- RunAzimuth(TP522, reference = "heartref")
Prediction_TP584 <- RunAzimuth(TP584, reference = "heartref")


## Change names on the rows of the Predictions seurat objects to nuclei ID + Animal ID#
rownames_TP31 <- rownames(Prediction_TP31@meta.data)
rownames_TP34 <- rownames(Prediction_TP34@meta.data)
rownames_TP41 <- rownames(Prediction_TP41@meta.data)
rownames_TP24 <- rownames(Prediction_TP24@meta.data)
rownames_TP522 <- rownames(Prediction_TP522@meta.data)
rownames_TP584 <- rownames(Prediction_TP584@meta.data)

new_rownames_TP31 <- gsub("-1$", "-TP31", rownames_TP31)
new_rownames_TP34 <- gsub("-1$", "-TP34", rownames_TP34)
new_rownames_TP41 <- gsub("-1$", "-TP41", rownames_TP41)
new_rownames_TP24 <- gsub("-1$", "-TP24", rownames_TP24)
new_rownames_TP522 <- gsub("-1$", "-TP522", rownames_TP522)
new_rownames_TP584 <- gsub("-1$", "-TP584", rownames_TP584)

colnames(Prediction_TP31) <- new_rownames_TP31
colnames(Prediction_TP34) <- new_rownames_TP34
colnames(Prediction_TP41) <- new_rownames_TP41
colnames(Prediction_TP24) <- new_rownames_TP24
colnames(Prediction_TP522) <- new_rownames_TP522
colnames(Prediction_TP584) <- new_rownames_TP584

## Merge all 6 animals into one Seurat object ##

mouse_1w <- merge(Prediction_TP31, 
                  y=c(Prediction_TP34, 
                      Prediction_TP41, 
                      Prediction_TP24, 
                      Prediction_TP522, 
                      Prediction_TP584), 
                  project = "mouse_1w_full_merge_project")
mouse_1w <- NormalizeData(mouse_1w, 
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000)
mouse_1w <- FindVariableFeatures(mouse_1w)
mouse_1w <- ScaleData(mouse_1w)
mouse_1w <- RunPCA(mouse_1w, npcs = 20)
mouse_1w <- RunUMAP(mouse_1w, dims = 1:20)
mouse_1w <- IntegrateLayers(object = mouse_1w, 
                                      method = RPCAIntegration, 
                                      orig.reduction = "pca", 
                                      new.reduction = "integrated.rpca",
                            verbose = FALSE)
mouse_1w <- JoinLayers(mouse_1w)


# Save Seurat object
SaveH5Seurat(mouse_1w, 
             filename = "mouse_1w.h5seurat",
             overwrite = TRUE)


## Filter for ventricular cardiomyocytes?
mouse_1w_vcm <- subset(mouse_1w, 
                          predicted.celltype.l2.score >= 0.7 & 
                            mapping.score >= 0.7 & 
                            predicted.celltype.l1.score >= 0.7 & 
                            mapping.score >= 0.7 &
                            predicted.celltype.l2 == 'Ventricular Cardiomycoyte')
head(mouse_1w_vcm, 5) #control post: inspect the top 5 rows in the dataset
tail(mouse_1w_vcm, 5)


# Assign the column "orig.ident" as the identity columns for the Seurat object.
Idents(mouse_1w_vcm) <- mouse_1w_vcm@meta.data$orig.ident

# Find differentially expressed genes
mouse_1w_vcm_markers <- FindMarkers(mouse_1w_vcm, ident.1 = "1 Week - SHAM", ident.2 = "1 Week - AB")
    
mouse_1w_vcm_markers_all_genes <- mouse_1w_vcm_markers %>%
    tibble::rownames_to_column(var = "gene") %>%
    mutate(log_pval = -log10(p_val_adj))
write.xlsx(mouse_1w_vcm_markers_all_genes, file = "mouse_1w_vcm_markers_all_genes.xlsx")

# Filter for significant genes
mouse_1w_vcm_markers_significant_genes <- mouse_1w_vcm_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
write.xlsx(mouse_1w_vcm_markers_significant_genes, file = "mouse_1w_vcm_markers_significant_genes.xlsx") # Save data frame as excel