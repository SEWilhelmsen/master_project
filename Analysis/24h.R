### RNA sequencing data SHAM vs. AB - Preprocessing
## Timepoint: 24 hours
## SHAM: 'TP17-2'
## AB: 'TP14-2'
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



## Load the data - 1 sham and 1 AB mice - 24 hours after ORAB operation

counts <- Seurat::Read10X_h5("W:/CardioTarget/Multiome seq data/TP17-2/outs/filtered_feature_bc_matrix.h5")
TP172 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "24 Hours - SHAM"
)

counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/TP14-2/outs/filtered_feature_bc_matrix.h5")
TP142 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "24 Hours - AB"
)

## Annotation of nuclei - Assignment to each nuclei a score on what type of cell type the nucleus comes from

options(future.globals.maxSize = 5 * 1024^3) # If error on the next lines, change max size to 2 GB. 
getOption("future.globals.maxSize")

Prediction_TP172 <- RunAzimuth(TP172, reference = "heartref")
Prediction_TP142 <- RunAzimuth(TP142, reference = "heartref")


## Change names on the rows of the Predictions seurat objects to nuclei ID + Animal ID#
rownames_TP172 <- rownames(Prediction_TP172@meta.data)
rownames_TP142 <- rownames(Prediction_TP142@meta.data)

new_rownames_TP172 <- gsub("-1$", "-TP172", rownames_TP172)
new_rownames_TP142 <- gsub("-1$", "-TP142", rownames_TP142)

colnames(Prediction_TP172) <- new_rownames_TP172
colnames(Prediction_TP142) <- new_rownames_TP142

# Run PCA on each condition

mouse_24h_sham_normalize <- NormalizeData(Prediction_TP172, normalization.method = "LogNormalize", scale.factor = 10000)
mouse_24h_sham_find_features <- FindVariableFeatures(mouse_24h_sham_normalize)
mouse_24h_sham_scale <- ScaleData(mouse_24h_sham_find_features)
mouse_24h_sham_pca <- RunPCA(mouse_24h_sham_scale, npcs = 20)

# Plot of clusters is C8glei. Should the data be more specific? 
DimPlot(mouse_24h_sham_pca, reduction = "pca")

## Merge all 4 animals into one Seurat object

mouse_24h <- merge(Prediction_TP172, y=c(Prediction_TP142), project = "Full_Merge_Project_24h")
mouse_24h_normalize <- NormalizeData(mouse_24h, normalization.method = "LogNormalize", scale.factor = 10000)
mouse_24h_find_features <- FindVariableFeatures(mouse_24h_normalize)
mouse_24h_scale <- ScaleData(mouse_24h_find_features)
mouse_24h_pca <- RunPCA(mouse_24h_scale, npcs = 20)

# Plot of clusters is C8glei. Should the data be more specific? 
DimPlot(mouse_24h_pca, reduction = "pca")

mouse_24h_umap <- RunUMAP(mouse_24h_pca, dims = 1:20)
mouse_24h_integrate <- IntegrateLayers(object = mouse_24h_umap, 
                             method = RPCAIntegration, 
                             orig.reduction = "pca", 
                             new.reduction = "integrated.rpca",
                            verbose = FALSE)
mouse_24h_join_layers <- JoinLayers(mouse_24h_integrate)


# Save Seurat object
SaveH5Seurat(mouse_24h_join_layers, overwrite = TRUE)

mouse_24h_join_layers_assay <- GetAssayData(object = mouse_24h_join_layers, assay = "RNA", layer = "data")
write.csv(mouse_24h_join_layers_assay, "mouse_24h_join_layers.csv") #this takes a looong time

## Start here to directly load the Seurat object
# Load Seurat object
hfile <- Connect("Full_Merge_Project_24h.h5Seurat")
hfile


## Filter for ventricular cardiomyocytes?
mouse_24h_vcm <- subset(mouse_24h_join_layers, 
                          predicted.celltype.l2.score >= 0.7 & 
                            mapping.score >= 0.7 & 
                            predicted.celltype.l1.score >= 0.7 & 
                            mapping.score >= 0.7 &
                            predicted.celltype.l2 == 'Ventricular Cardiomycoyte')
head(mouse_24h_vcm, 5) #control post: inspect the top 5 rows in the dataset
tail(mouse_24h_vcm, 5)

SaveH5Seurat(mouse_24h_vcm, "C:/Users/siljeew/snRNAseq/24h_sham_AB/mouse_24h_vcm.h5seurat", overwrite = TRUE)

## Volcano plot 
# Assign the column "orig.ident" as the identity columns for the Seurat object.
Idents(mouse_24h_vcm) <- mouse_24h_vcm@meta.data$orig.ident

# Find differentially expressed genes
markers <- FindMarkers(mouse_24h_vcm, ident.1 = "24 Hours - SHAM", ident.2 = "24 Hours - AB")
    
markers <- markers %>%
    tibble::rownames_to_column(var = "gene") %>%
    mutate(log_pval = -log10(p_val_adj))
    
top_genes_200 <- markers %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    head(200)

top_genes_VCM_24h <- markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC))) 

write_xlsx(top_genes_VCM_24h, "C:/Users/siljeew/snRNAseq/24h_sham_AB/top_genes_VCM_24h.xlsx") # Save data frame as excel


# 0: AB = SHAM, +: mer i AB enn SHAM, -: mindre i AB enn SHAM
ggplot(markers, aes(x = avg_log2FC, y = log_pval)) +
  geom_point(aes(color = p_val_adj < 0.05)) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot VCM - 24 Hours",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "none") +
  geom_text_repel(data = top_genes_200, aes(label = gene), size = 3) +
  ylim(0, 50)



## GO gene ontology 

markers$RowNames <- rownames(markers)
markers <- FindAllMarkers(mouse_24h_vcm)
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

barplot(go_results_VCM,showCategory = 20, title = "20 top biological processes - 24 Hours")
dotplot(go_results_VCM, showCategory = 20, title = "20 top biological processes - 24 Hours")

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
p1 <- DimPlot(Mouse_24h, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p1

p2 <- DimPlot(Mouse_24h_VCM, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p2

head(Mouse_24h_VCM)