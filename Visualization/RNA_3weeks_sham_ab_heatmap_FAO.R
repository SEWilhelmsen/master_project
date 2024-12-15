## Heatmap 
# 1. Filter genes for specific metabolic process
# 2. Subset conditions 
# 3. Find mean of samples
# 4. Create heatmap

library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)


##### 1. Filter genes for specific metabolic process

# Keep only significant genes
sign_genes <- markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))

head(sign_genes) # Control post

#### 2. Define GO terms and find related genes 

# Define GO terms
GO_ID <- list(
  FAO = "GO:0019395"
)

# Function to retrieve symbols for a GO term
get_genes_for_GO <- function(GO_ID) {
  entrez_ids <- AnnotationDbi::select(org.Hs.eg.db, 
                                      keys = GO_ID, 
                                      keytype = "GOALL", 
                                      columns = "ENTREZID")[, "ENTREZID"]
  gene_symbols <- AnnotationDbi::select(org.Hs.eg.db, 
                                        keys = entrez_ids, 
                                        keytype = "ENTREZID", 
                                        columns = "SYMBOL")[, "SYMBOL"]
  gene_symbols
}


GO_genes <- lapply(GO_ID, get_genes_for_GO) # Retrieve genes for each GO term

#### 3. Filter GO genes in the dataset
# Normalize gene symbols to ensure they match case and format
normalize_symbols <- function(genes) {
  toupper(genes)
}
normalized_sign_genes <- normalize_symbols(sign_genes$gene)
normalized_GO_genes <- lapply(GO_genes, normalize_symbols)

all_GO_genes <- unique(unlist(normalized_GO_genes)) # Unlist GO_genes and create a unique list of GO related genes
common_genes <- intersect(normalized_sign_genes, all_GO_genes) # Find the intersection between sign_genes and GO_genes
sign_genes_subset <- sign_genes[sign_genes$gene %in% common_genes, ] # Subset sign_genes based on common_genes


rownames(sign_genes_subset) <- sign_genes_subset$gene # Convert the data frame to a matrix
sign_genes_numeric <- sign_genes_subset[, sapply(sign_genes_subset, is.numeric)] # Remove the gene column to keep only numerical data

# Convert to a matrix
sign_genes_matrix <- as.matrix(sign_genes_numeric)
sign_genes_matrix <- sign_genes_matrix[, "avg_log2FC", drop = FALSE]
colnames(sign_genes_matrix)[1] <- "3 weeks" # Change the name of column 1 to "3 weeks"


# Reorder the matrix based on GO terms
gene_processes <- data.frame(
  gene = unlist(normalized_GO_genes),
  process = rep(names(GO_ID), times = sapply(normalized_GO_genes, length))
)
# Merge with sign_genes to ensure all common genes are considered and create a sorted table based on GO processes
merged_data <- merge(data.frame(gene = rownames(sign_genes_matrix)), gene_processes, by = "gene", all.x = TRUE)

merged_data$process[is.na(merged_data$process)] <- "Unknown" # Fill process NA values with "Unknown" or another category
merged_data <- merged_data %>% arrange(process) # Sort data by process
merged_data <- merged_data %>% # Remove possible duplicates in merged data
  arrange(process, gene) %>%
  distinct(gene, .keep_all = TRUE)


sorted_genes <- merged_data$gene  # sorted_genes should be a character vector
sign_genes_matrix <- sign_genes_matrix[sorted_genes, , drop = FALSE] # Reorder the matrix rows based on sorted_genes

# Verify dimensions
cat("Number of rows in sign_genes_matrix: ", nrow(sign_genes_matrix), "\n")
cat("Number of genes in row annotations: ", length(merged_data$process), "\n")


row_anno <- rowAnnotation(
  Process = merged_data$process,
  col = list(Process = c("FAO" = "aquamarine3")),
  show_annotation_name = TRUE
)

### 4. Create a heatmap of avg_log2FC
heatmap <- Heatmap(sign_genes_matrix,
                   name = "Expression",
                   cluster_rows = TRUE, 
                   show_row_names = TRUE, 
                   row_names_gp = gpar(fontsize = 8),
                   row_names_side = "left",
                   row_title = "Genes",
                   show_row_dend = FALSE,
                   left_annotation = row_anno,
                   cluster_columns = TRUE,  
                   show_column_names = TRUE,
                   column_title = "Expression AB vs SHAM",
                   show_column_dend = FALSE,
                   heatmap_legend_param = list(
                     title = "avg_log2FC",
                     legend_direction = "vertical",
                     legend_height = unit(4, "cm")
                   ),
                   gap = unit(1, "mm"),
                   heatmap_height = unit(0.4, "cm")*nrow(sign_genes_matrix),
                   col = colorRamp2(c(min(-1), 0, max(1)),
                                    c("navy", "white", "red"))
)

draw(heatmap)