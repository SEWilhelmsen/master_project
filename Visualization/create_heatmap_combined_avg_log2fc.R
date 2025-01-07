## Create heatmap from combined_avg_log2fc
# Silje Wilhelmsen

# Function to create the gene process mapping
create_gene_process_mapping <- function(genes_of_interest_for_heatmap_vector) {
  gene_process_mapping <- list()
  
  # Ensure correct indexing on the list
  for (process in names(genes_of_interest_for_heatmap_vector)) {
    genes <- genes_of_interest_for_heatmap_vector[[process]]  # ensure this accesses correct index
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


# Prepare heatmap matrix should consider global variable dependency
prepare_heatmap_matrix <- function(data, genes_of_interest_for_heatmap_vector) {
  # Ensure both lists are uppercase
  genes_of_interest_for_heatmap_vector <- toupper(genes_of_interest_for_heatmap_vector)
  
  # Subset data
  subset <- data[rownames(data) %in% genes_of_interest_for_heatmap_vector, ]
  
  print("Subset after filtering:")
  print(subset)
  
  # Add process information
  subset$go_processes <- sapply(rownames(subset), function(gene) {
    gene <- toupper(gene)
    processes <- gene_process_mapping[[gene]]
    if (!is.null(processes)) {
      paste(processes, collapse = ",")
    } else {
      "None"
    }
  })
  
  as.matrix(subset[, 1:6])
}


# Function to create heatmap
create_heatmap <- function(heatmap_matrix_numeric, go_association_vector) {
  if (!is.matrix(heatmap_matrix_numeric)) {
    stop("Input data must be a matrix.")
  }
  
  heatmap_matrix_numeric[is.na(heatmap_matrix_numeric)] <- 0
  heatmap_matrix_numeric[!is.finite(heatmap_matrix_numeric)] <- 0
  
  if (nrow(heatmap_matrix_numeric) != length(go_association_vector)) {
    stop("Mismatch between number of matrix rows and association vector length.")
  }
  
  # Define color function
  col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  # Create text annotation with non-empty vector
  row_annotation <- rowAnnotation(
    text = anno_text(go_association_vector, gp = gpar(fontsize = 8)),
    show_annotation_name = FALSE
  )
  
  # Create the heatmap
  heatmap <- Heatmap(
    heatmap_matrix_numeric, 
    name = "log2 fold change",
    col = col_fun,
    na_col = "grey",
    cluster_rows = TRUE, 
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 10),
    row_gap = unit(3, "mm"),
    cluster_columns = FALSE,  
    show_column_names = TRUE,
    column_names_rot = 45,
    column_names_centered = TRUE,
    column_title = "Differential Gene Expression AB vs SHAM",
    show_column_dend = FALSE,
    width = unit(10, "cm"),
    heatmap_legend_param = list(
      title = "avg_log2fc",
      legend_direction = "vertical",
      legend_height = unit(4, "cm")
    ),
    right_annotation = row_annotation,  # Use right_annotation for adding text
    gap = unit(1, "mm")
  )
  
  # Create custom legend for NA
  na_legend <- Legend(
    labels = "NA",
    title = "",
    legend_gp = gpar(fill = "grey"),
    direction = "vertical"
  )
  
  # Draw the heatmap and stack legends
  draw(heatmap, annotation_legend_list = list(na_legend), merge_legend = TRUE)
}