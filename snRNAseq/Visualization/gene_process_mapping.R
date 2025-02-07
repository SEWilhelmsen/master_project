# Gene mapping for plots
# Silje Wilhelmsen

# Function to create the gene process mapping
create_gene_process_mapping <- function(genes_of_interest_for_heatmap_vector) {
  gene_process_mapping <- list()
  
  # Ensure correct indexing on the list
  for (process in names(genes_of_interest_for_heatmap_vector)) {
    genes <- genes_of_interest_for_heatmap_vector[[process]]
    for (gene in genes) {
      # Check symbol normalization for consistency
      gene <- toupper(gene)  # Normalize (if required) to match rownames format
      if (!is.null(gene_process_mapping[[gene]])) {
        gene_process_mapping[[gene]] <- c(gene_process_mapping[[gene]], process)
      } else {
        gene_process_mapping[[gene]] <- process
      }
    }
  }
  
  gene_process_mapping
}

