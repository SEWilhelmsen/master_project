# Functions to create network plot 
# Silje Wilhelmsen

# Define function to get genes of interest
get_genes_of_interest_vector <- function(go_processes_of_interest, define_go_process, get_genes_for_go_process) {
  genes_of_interest_vector <- get_genes_of_interest(go_processes_of_interest, define_go_process, get_genes_for_go_process)
  return(unlist(genes_of_interest_vector, use.names = FALSE))
}

# Function to load and preprocess data
load_and_preprocess_data <- function(file_path, go_processes_of_interest, define_go_process, get_genes_for_go_process) {
  mouse_6h_vcm <- readRDS(file_path)
  genes_of_interest_vector <- get_genes_of_interest_vector(go_processes_of_interest, define_go_process, get_genes_for_go_process)
  
  expression_matrix <- GetAssayData(mouse_6h_vcm, layer = "data")
  dense_expression_matrix <- as.matrix(expression_matrix)
  subset_expression_matrix <- dense_expression_matrix[rownames(dense_expression_matrix) %in% genes_of_interest_vector, ]
  
  return(subset_expression_matrix)
}

# Function to compute correlation matrix
compute_correlation_matrix <- function(subset_expr_matrix) {
  zero_variance_genes <- rownames(subset_expr_matrix)[apply(subset_expr_matrix, 1, sd) == 0]
  if (length(zero_variance_genes) > 0) {
    subset_expr_matrix <- subset_expr_matrix[!rownames(subset_expr_matrix) %in% zero_variance_genes, ]
  }
  
  if (nrow(subset_expr_matrix) > 0) {
    correlation_matrix <- cor(t(subset_expr_matrix), use = "pairwise.complete.obs")
    return(correlation_matrix)
  } else {
    stop("No data left after filtering zero-variance genes.")
  }
}

# Create adjacency matrix from correlation
create_adjacency_matrix <- function(correlation_matrix, genes, threshold) {
  n <- length(genes)
  adjacency_matrix <- matrix(0, nrow = n, ncol = n, dimnames = list(genes, genes))
  
  if (n > 1) {
    for (i in 1:n) {
      for (j in i:n) {
        if (!is.na(correlation_matrix[i, j]) && abs(correlation_matrix[i, j]) > threshold && i != j) {
          adjacency_matrix[i, j] <- 1
          adjacency_matrix[j, i] <- 1  # Undirected
        }
      }
    }
  } else {
    warning("Correlation matrix is too small for adjacency computation.")
  }
  return(adjacency_matrix)
}

# Function to retrieve GO terms for a given list of gene symbols
get_go_terms_for_genes <- function(gene_symbols) {
  entrez_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbols, keytype = "SYMBOL", columns = "ENTREZID")[, "ENTREZID"]
  gene_go_pairs <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez_ids, keytype = "ENTREZID", columns = c("SYMBOL", "GOALL"))
  return(gene_go_pairs)
}