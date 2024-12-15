### RNA sequencing data SHAM vs. AB - Preprocessing
## Timepoint: 1 Week
## SHAM: 'TP3-1', 'TP3-4', 'TP4-1'
## AB: 'TP2-4', 'TP52-2', 'TP58-4'


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
library(devtools)

# AvailableData()

## Load the data - 3 sham and 3 AB mice - 1 week after ORAB operation
## Error in local. Crsparse ... not found -> Try leave the laptop off for a night, suddenly it work
## Can also try to switch computer, that might work too

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

# SaveH5Seurat(TP31, file = "C:/Users/siljeew/snRNAseq/1week_sham_AB/TP31.h5seurat", overwrite = TRUE)
# SaveH5Seurat(TP34, file = "C:/Users/siljeew/snRNAseq/1week_sham_AB/TP34.h5seurat", overwrite = TRUE)
# SaveH5Seurat(TP41, file = "C:/Users/siljeew/snRNAseq/1week_sham_AB/TP41.h5seurat", overwrite = TRUE)
# SaveH5Seurat(TP24, file = "C:/Users/siljeew/snRNAseq/1week_sham_AB/TP24.h5seurat", overwrite = TRUE)
# SaveH5Seurat(TP522, file = "C:/Users/siljeew/snRNAseq/1week_sham_AB/TP522.h5seurat", overwrite = TRUE)
# SaveH5Seurat(TP584, file = "C:/Users/siljeew/snRNAseq/1week_sham_AB/TP584.h5seurat", overwrite = TRUE)


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
                  project = "mouse_1w_Full_Merge_Project")
mouse_1w_normalize <- NormalizeData(mouse_1w, 
                                    normalization.method = "LogNormalize", 
                                    scale.factor = 10000)
mouse_1w_find_features <- FindVariableFeatures(mouse_1w_normalize)
mouse_1w_scale <- ScaleData(mouse_1w_find_features)
mouse_1w_pca <- RunPCA(mouse_1w_scale, npcs = 20)
mouse_1w_umap <- RunUMAP(mouse_1w_pca, dims = 1:20)
mouse_1w_integrate <- IntegrateLayers(object = mouse_1w_umap, 
                                      method = RPCAIntegration, 
                                      orig.reduction = "pca", 
                                      new.reduction = "integrated.rpca",
                            verbose = FALSE)
mouse_1w_join_layers <- JoinLayers(mouse_1w_integrate)


# Save Seurat object
SaveH5Seurat(mouse_1w_join_layers, 
             filename = "mouse_1w_join_layers.h5seurat",
             overwrite = TRUE)


## Filter for ventricular cardiomyocytes?
mouse_1w_vcm <- subset(mouse_1w_join_layers, 
                          predicted.celltype.l2.score >= 0.7 & 
                            mapping.score >= 0.7 & 
                            predicted.celltype.l1.score >= 0.7 & 
                            mapping.score >= 0.7 &
                            predicted.celltype.l2 == 'Ventricular Cardiomycoyte')
head(mouse_1w_vcm, 5) #control post: inspect the top 5 rows in the dataset
tail(mouse_1w_vcm, 5)


## Volcano plot 
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

# Filter for top 200 significant genes
mouse_1w_vcm_markers_significant_genes_top_200 <- mouse_1w_vcm_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(200)
write.xlsx(mouse_1w_vcm_markers_significant_genes_top_200, file = "mouse_1w_vcm_markers_significant_genes_top_200.xlsx") # Save data frame as excel



########### Dont know what is underneath here
# 0: AB = SHAM, +: mer i AB enn SHAM, -: mindre i AB enn SHAM
ggplot(markers, aes(x = avg_log2FC, y = log_pval)) +
  geom_point(aes(color = p_val_adj < 0.05)) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot VCM - 1 week",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "none") +
  geom_text_repel(data = top_genes_200, aes(label = gene), size = 3) +
  ylim(0, 50)




## GO gene ontology 

markers$RowNames <- rownames(markers)
markers <- FindAllMarkers(mouse_1w_VCM)
significant_genes <- markers[markers$p_val_adj < 0.05, "gene"]

gene_list <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform GO analysis
go_results_VCM <- enrichGO(gene = gene_list$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",  # Biological Process
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

barplot(go_results_VCM,showCategory = 20, title = "Top 20 processes of top expressed genes?- 1 Week")
dotplot(go_results_VCM, showCategory = 20, title = "Top 20 processes of top expressed genes?- 1 Week")

go_results_df <- as.data.frame(go_results_VCM)
go_results_df_sorted <- go_results_df[order(go_results_df$p.adjust), ]
head(go_results_df_sorted) # This gave the same result as go_results_df. Is this automatically done in the previous code? 

nrow(go_results_df) #Inspect the number of rows/observations in the data frame
ncol(go_results_df) #Inspect the number of columns in the data frame


## Search: "lipid"
lipid_modification_results <- go_results_df[grep("lipid", go_results_df$Description, ignore.case = TRUE), ]

# Sort by p.adjust in descending order
lipid_modification_results <- lipid_modification_results %>%
  mutate(Description = reorder(Description, -p.adjust))
# Create a bar plot with sorted metabolic processes
bar_plot <- ggplot(lipid_modification_results, +
                   aes(x = Description, y = p.adjust)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(x = "Metabolism Term", y = "p.adjust", title = "Significant Lipid Modification Results within VCM")
# Display the bar plot
print(bar_plot)
# View the filtered results
print(lipid_modification_results)


## Search: "fatty acid"
FA_results <- go_results_df[grep("fatty acid", go_results_df$Description, ignore.case = TRUE), ]
FA_results <- FA_results %>%
  mutate(Description = reorder(Description, -p.adjust))
bar_plot <- ggplot(FA_results, aes(x = Description, y = p.adjust)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(x = "Metabolism Term", y = "p.adjust", title = "Significant Fatty Acid Results within VCM")
print(bar_plot)
print(FA_results)



## Search: "metabo"
metabo_results <- go_results_df[grep("metabo", go_results_df$Description, ignore.case = TRUE), ]
# Sort the metabolic processes by p.adjust in descending order
metabo_results_sorted <- metabo_results %>%
  mutate(Description = reorder(Description, -p.adjust))
# Create a bar plot with sorted metabolic processes
bar_plot <- ggplot(metabo_results_sorted, aes(x = Description, y = p.adjust)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(x = "Metabolism Term", y = "p.adjust", title = "Significant Metabolism Results within VCM")
# Display the bar plot
print(bar_plot)
print(metabo_results)

## Search: "carb"
carb_results <- go_results_df[grep("carb", go_results_df$Description, ignore.case = TRUE), ]
carb_results <- carb_results %>%
  mutate(Description = reorder(Description, -p.adjust))
bar_plot <- ggplot(carb_results, aes(x = Description, y = p.adjust)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(x = "Metabolism Term", y = "p.adjust", title = "Significant Carb Results within VCM")
print(bar_plot)
print(carb_results)

## Search: "beta"
beta_results <- go_results_df[grep("beta", go_results_df$Description, ignore.case = TRUE), ]
beta_results <- beta_results %>%
  mutate(Description = reorder(Description, -p.adjust))
bar_plot <- ggplot(beta_results, aes(x = Description, y = p.adjust)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(x = "Metabolism Term", y = "p.adjust", title = "Significant beta Results within VCM")
print(bar_plot)
print(beta_results)

## Search: "keto"
keto_results <- go_results_df[grep("keto", go_results_df$Description, ignore.case = TRUE), ]
keto_results <- keto_results %>%
  mutate(Description = reorder(Description, -p.adjust))
bar_plot <- ggplot(keto_results, aes(x = Description, y = p.adjust)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(x = "Metabolism Term", y = "p.adjust", title = "Significant keto Results within VCM")
print(bar_plot)
print(keto_results)

## Visualize the output of Azimuth
p1 <- DimPlot(mouse_1w, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p1

p2 <- DimPlot(mouse_1w_VCM, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p2

head(mouse_1w_VCM)