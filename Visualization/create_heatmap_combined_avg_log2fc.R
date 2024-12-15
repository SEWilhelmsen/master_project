## Create heatmap from combined_avg_log2fcHEATMAP
# Silje Wilhelmsen
# 1. Filter significant genes in the dataset
# 2. Define GO genes
# 3. Filter GO genes in the dataset
# 4. Create a heatmap



library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(ComplexHeatmap) #https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#row-and_column_orders
library(circlize)
library(RColorBrewer)
library(readr)
library(tibble)

## Load files
# Change the name of the first column manually to "gene" prior to opening the .csv here
common_genes_with_avg_log2fc <- read_csv("common_genes_with_avg_log2fc.csv")

### Define GO terms and find related genes 
source("define_go_terms_and_fins_genes.R")


# Subset the data frame to include only the genes in all_go_genes
filtered_common_genes <- common_genes_with_avg_log2FC %>%
  filter(gene %in% all_go_genes)

# Create an identifier column for each GO category and determine their corresponding genes
filtered_common_genes <- filtered_common_genes %>%
  mutate(go_category = case_when(
    gene %in% normalized_go_genes$glycolysis ~ 'Glycolysis',
    gene %in% normalized_go_genes$fao ~ 'FAO',
    gene %in% normalized_go_genes$tca ~ 'TCA'
  ))

print(head(filtered_common_genes)) # Check the first few rows of the filtered data frame


# Prepare the data for heatmap creation
gene_matrix <- filtered_common_genes %>%
  column_to_rownames(var = "gene")


# Convert the subset of data to a matrix
gene_matrix_ordered <- as.matrix(gene_matrix[, c("avg_log2FC_24h", "avg_log2FC_3d", "avg_log2FC_1w", "avg_log2FC_3w")])

# Print the first few rows of the ordered gene matrix to confirm the column order
print(colnames(gene_matrix_ordered))
print(head(gene_matrix_ordered))

# Create a heatmap sorted by GO category
Heatmap(gene_matrix_ordered, 
        name = "log2fc", 
        cluster_rows = FALSE, 
        show_row_names = TRUE,
        row_names_side = "left",
        row_split = filtered_common_genes$go_category,
        row_names_gp = gpar(fontsize = 12),
        row_gap = unit(2, "mm"),
        cluster_columns = FALSE,  
        show_column_names = TRUE,
        column_names_rot = 30,                          # Rotate column names 30*
        column_title = "Expression AB vs SHAM",
        show_column_dend = FALSE,
        heatmap_legend_param = list(
          title = "avg_log2FC",
          legend_direction = "vertical",
          legend_height = unit(4, "cm")
        ),
        gap = unit(2, "mm"),
        heatmap_height = unit(5, "cm")*nrow(gene_matrix_ordered),
        col = colorRamp2(c(min(-1), 0, max(1)),
                         c("navy", "white", "red4"))
)
