## Find significant genes 
## Output: table of significant genes list??
# Silje Wilhelmsen

# Keep only significant genes
### Define Function to Process Markers for Each Time Point ###
process_markers <- function(markers) {
  sign_genes <- markers %>%
    filter(p_val_adj < 0.05) %>%
    arrange(desc(abs(avg_log2FC)))
  
  normalized_sign_genes <- toupper(sign_genes$gene)
  
  list(sign_genes = sign_genes, normalized_sign_genes = normalized_sign_genes)
}

# Process markers for each time point
markers_3w <- process_markers(markers_3w)
markers_1w <- process_markers(markers_1w)
markers_3d <- process_markers(markers_3d)