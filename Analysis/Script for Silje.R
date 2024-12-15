### Script for Silje
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

#install.packages("heartref.SeuratData", repos = "http://seurat.nygenome.org", type = "source")
library(heartref.SeuratData)


# AvailableData()

## Load the data - 2 sham and 2 AB mice - 3 weeks after ORAB operation

counts <- Seurat::Read10X_h5("W:/CardioTarget/Multiome seq data/GCF0742_GCF0795_CellRanger/Sample2-TP511B/outs/filtered_feature_bc_matrix.h5")
TP511 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "3 Weeks - AB"
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

counts <- Read10X_h5("W:/CardioTarget/Multiome seq data/GCF0742_GCF0795_CellRanger/Sample4-TP62A/outs/filtered_feature_bc_matrix.h5")
TP62 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA", project = "3 Weeks - AB"
)

## Annotation of nuclei - Assignment to each nuclei a score on what type of cell type the nucleus comes from

Prediction_TP62 <- RunAzimuth(TP62, reference = "heartref")
Prediction_TP511 <- RunAzimuth(TP511, reference = "heartref")
Prediction_TP52 <- RunAzimuth(TP52, reference = "heartref")
Prediction_TP75 <- RunAzimuth(TP75, reference = "heartref")


## Change names on the rows of the Predictions seurat objects to nuclei ID + Animal ID#
rownames_TP511 <- rownames(Prediction_TP511@meta.data)
rownames_TP62 <- rownames(Prediction_TP62@meta.data)
rownames_TP75 <- rownames(Prediction_TP75@meta.data)
rownames_TP52 <- rownames(Prediction_TP52@meta.data)

new_rownames_TP511 <- gsub("-1$", "-TP511", rownames_TP511)
new_rownames_TP62 <- gsub("-1$", "-TP62", rownames_TP62)
new_rownames_TP75 <- gsub("-1$", "-TP75", rownames_TP75)
new_rownames_TP52 <- gsub("-1$", "-TP52", rownames_TP52)

colnames(Prediction_TP511)<-new_rownames_TP511
colnames(Prediction_TP62)<-new_rownames_TP62
colnames(Prediction_TP75)<-new_rownames_TP75
colnames(Prediction_TP52)<-new_rownames_TP52

## Merge all 4 animals into one Seurat object 

Mouse_3w <- merge(Prediction_TP62, y=c(Prediction_TP511, Prediction_TP52, Prediction_TP75), project = "Full Merge Project")
Mouse_3w <- NormalizeData(Mouse_3w, normalization.method = "LogNormalize", scale.factor = 10000)
Mouse_3w <- FindVariableFeatures(Mouse_3w)
Mouse_3w <- ScaleData(Mouse_3w)
Mouse_3w <- RunPCA(Mouse_3w, npcs = 20)
Mouse_3w <- RunUMAP(Mouse_3w, dims = 1:20)
Mouse_3w <- IntegrateLayers(object = Mouse_3w, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                            verbose = FALSE)
Mouse_3w <- JoinLayers(Mouse_3w)


# Save Seurat object
SaveH5Seurat(Mouse_3w, overwrite = TRUE)

## Start here to directly load the Seurat object
# Load Seurat object
hfile <- Connect("Full Merge Project.h5Seurat")
hfile


## QC
Mouse_3w_filter <- subset(Mouse_3w, subset = predicted.celltype.l2.score >= 0.5 & mapping.score >= 0.5 & predicted.celltype.l1.score >= 0.5 & mapping.score >= 0.5)


## Volcano plot 
# Assign the column "orig.ident" as the identity columns for the Seurat object.
Idents(Mouse_3w) <- Mouse_3w@meta.data$orig.ident

# Find differentially expressed genes
markers <- FindMarkers(Mouse_3w, ident.1 = "3 Weeks - AB", ident.2 = "3 Weeks - SHAM")
    
markers <- markers %>%
    tibble::rownames_to_column(var = "gene") %>%
    mutate(log_pval = -log10(p_val_adj))
    
top_genes <- markers %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    head(200)


# 0: AB = SHAM, +: mer i AB enn SHAM, -: mindre i AB enn SHAM
ggplot(markers, aes(x = avg_log2FC, y = log_pval)) +
  geom_point(aes(color = p_val_adj < 0.05)) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "none") +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3) +
  ylim(0, 50)


## GO gene ontology 

markers$RowNames <- rownames(markers)
markers <- FindAllMarkers(Mouse_3w)
significant_genes <- markers[markers$p_val_adj < 0.05, "gene"]

gene_list <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform GO analysis
go_results <- enrichGO(gene = gene_list$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",  # Biological Process
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

barplot(go_results, showCategory = 20)
dotplot(go_results, showCategory = 20)


go_results_df <- as.data.frame(go_results)



# ## GO analysis in "Ventricular Cardiomyocyte"
# 
# # Obtain the "Ventricular Cardiomyocyte" subgroup based on predicted cell types
# cardiomyocyte_subgroup <- subset(Mouse_3w, predicted.celltype.l2 == "Ventricular Cardiomycoyte")
# markers$RowNames <- rownames(markers) # Extract the row names of the markers
# markers <- FindAllMarkers(cardiomyocyte_subgroup) # Perform the FindAllMarkers function within the "Ventricular Cardiomyocyte" subgroup
# significant_genes <- markers[markers$p_val_adj < 0.05, "gene"] # Filter significant genes
# gene_list <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) # Convert gene symbols to ENTREZ IDs
# 
# # Perform GO analysis on "Ventricular Cardiomyocyte" genes
# go_results <- enrichGO(gene = gene_list$ENTREZID,
#                        OrgDb = org.Hs.eg.db,
#                        keyType = "ENTREZID",
#                        ont = "BP",  # Biological Process
#                        pAdjustMethod = "BH",
#                        pvalueCutoff = 0.05,
#                        qvalueCutoff = 0.05)
# 
# go_results_sorted <- go_results[order(go_results$qvalue, decreasing = TRUE), ] 
# # Transpose the data frame
# top_go_results_transposed <- t(top_go_results)
# 
# # Convert the values to numeric and plot a horizontal barplot for the top 20 enriched GO terms
# barplot(-log10(as.numeric(top_go_results_transposed[2, ])), horiz = TRUE, names.arg = top_go_results_transposed[1, ], las = 1, col = "skyblue", main = "Top 20 Enriched GO Terms", xlab = "-log10(qvalue)")



## Search: "lipid"
lipid_modification_results <- go_results_df[grep("lipid", go_results_df$Description, ignore.case = TRUE), ]

# Sort by p.adjust in descending order
lipid_modification_results <- lipid_modification_results %>%
  mutate(Description = reorder(Description, -p.adjust))
# Create a bar plot with sorted metabolic processes
bar_plot <- ggplot(lipid_modification_results_sorted, aes(x = Description, y = p.adjust)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(x = "Metabolism Term", y = "p.adjust", title = "Significant Lipid Modification Results")
# Display the bar plot
print(bar_plot)
# View the filtered results
print(lipid_modification_results)


## Search: "metabo"
metabo_results <- go_results_df[grep("metabo", go_results_df$Description, ignore.case = TRUE), ]
# Sort the metabolic processes by p.adjust in descending order
metabo_results_sorted <- metabo_results %>%
  mutate(Description = reorder(Description, -p.adjust))
# Create a bar plot with sorted metabolic processes
bar_plot <- ggplot(metabo_results_sorted, aes(x = Description, y = p.adjust)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(x = "Metabolism Term", y = "p.adjust", title = "Significant Metabolism Results")
# Display the bar plot
print(bar_plot)



# # Filter results for "Ventricular Cardiomyocyte" subgroup and "metabo" keyword
# metabo_results <- go_results_df %>%
#   filter(grepl("metabo", Description, ignore.case = TRUE) & predicted.celltype.l2 == "Ventricular Cardiomyocyte")
# 
# # Sort the metabolic processes by p.adjust in descending order
# metabo_results_sorted <- metabo_results %>%
#   mutate(Description = reorder(Description, -p.adjust))
# 
# # Create a bar plot with sorted metabolic processes
# bar_plot <- ggplot(metabo_results_sorted, aes(x = Description, y = p.adjust)) +
#   geom_col(fill = "red") +
#   coord_flip() +
#   labs(x = "Metabolism Term", y = "p.adjust", title = "Significant Metabolism Results for Ventricular Cardiomyocyte")
# 
# # Display the bar plot
# print(bar_plot)


# Table
table(Prediction_TP511@meta.data$predicted.celltype.l2)