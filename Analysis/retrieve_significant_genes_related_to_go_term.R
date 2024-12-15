## Retrieve significant genes related to gene ontology term
## Output: matrix of significant genes related to GO term
# Silje Wilhelmsen

# Normalize gene symbols
normalize_symbols <- function(genes) {
  toupper(genes)
}
normalized_GO_genes <- lapply(GO_genes, normalize_symbols)

all_GO_genes <- unique(unlist(normalized_GO_genes)) # Unlist GO_genes and create a unique list of GO related genes

# Function to filter significant genes based on GO terms
filter_GO_genes <- function(normalized_sign_genes, sign_genes) {
  common_genes <- intersect(normalized_sign_genes, all_GO_genes) # Find the intersection between sign_genes and GO_genes
  sign_genes %>% filter(toupper(gene) %in% common_genes) # Subset sign_genes based on common_genes
}

# Filter GO genes for each time point
sign_genes_3w_subset <- filter_GO_genes(markers_3w$normalized_sign_genes, markers_3w$sign_genes)
sign_genes_1w_subset <- filter_GO_genes(markers_1w$normalized_sign_genes, markers_1w$sign_genes)
sign_genes_3d_subset <- filter_GO_genes(markers_3d$normalized_sign_genes, markers_3d$sign_genes)

# Combine log fold change values into a single matrix
combined_logFC <- function(..., col_names) {
  combined <- Reduce(function(x, y) {
    merge(x, y, by = "gene", all = TRUE)
  }, list(...))
  
  colnames(combined) <- c("gene", col_names)
  combined
}

combined_sign_genes <- combined_logFC(
  data.frame(gene = rownames(sign_genes_3w_subset), `3 Weeks` = sign_genes_3w_subset$avg_log2FC),
  data.frame(gene = rownames(sign_genes_1w_subset), `1 Week` = sign_genes_1w_subset$avg_log2FC),
  data.frame(gene = rownames(sign_genes_3d_subset), `3 Days` = sign_genes_3d_subset$avg_log2FC),
  col_names = c("3 Weeks", "1 Week", "3 Days")
)

# Ensure row order is consistent
rownames(combined_sign_genes) <- combined_sign_genes$gene
combined_sign_genes <- combined_sign_genes[ , -1]

# Convert combined data frame to matrix
sign_genes_matrix <- as.matrix(combined_sign_genes)
sign_genes_matrix[is.na(sign_genes_matrix)] <- 0  # Replace NA with 0
