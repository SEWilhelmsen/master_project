## Create heatmap from combined_avg_log2fc
# Silje Wilhelmsen


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
create_heatmap <- function(heatmap_matrix_numeric, go_association_vector, process_colors) {
  if (!is.matrix(heatmap_matrix_numeric)) {
    stop("Input data must be a matrix.")
  }
  
  if (nrow(heatmap_matrix_numeric) != length(go_association_vector)) {
    stop("Mismatch between number of matrix rows and association vector length.")
  }
  
  # Define color function
  col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  # # Create text annotation with non-empty vector
  # row_annotation <- rowAnnotation(
  #   text = anno_text(go_association_vector, gp = gpar(fontsize = 8)),
  #   show_annotation_name = FALSE
  # )
  
  # Define annotations based on the process
  right_annotation <- rowAnnotation(
    process = anno_simple(go_association_vector, col = process_colors),
    show_annotation_name = FALSE
  )

  # Create the heatmap
  heatmap <- Heatmap(
    heatmap_matrix_numeric, 
    name = "log2 fold change",
    col = col_fun,
    na_col = "grey",
    cluster_rows = FALSE, 
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
      legend_height = unit(4, "cm"),
      at = c(-1, -0.5, 0, 0.5, 1),
      labels = c("-1", "Higher expression in SHAM",  "0", "Higher expression in ORAB", "1")
    ),
    right_annotation = right_annotation,
    gap = unit(1, "mm")
  )
  
  # Create custom legend for NA
  na_legend <- Legend(
    labels = "NA",
    title = "",
    legend_gp = gpar(fill = "grey"),
    direction = "vertical"
  )
  
  # Filter process_colors to only include those present in go_association_vector
  present_processes <- unique(go_association_vector)
  relevant_colors <- process_colors[names(process_colors) %in% present_processes]
  
  process_legend <- Legend(
    title = "Process",
    labels = names(relevant_colors),
    legend_gp = gpar(fill = relevant_colors),
    direction = "vertical"
  )
  
  # Draw the heatmap and stack all legends
  draw(heatmap, annotation_legend_list = list(na_legend, process_legend), merge_legend = TRUE)
}




# Function to save a heatmap using ComplexHeatmap and ensure filename reflects processes
save_heatmap_complex <- function(heatmap, output_dir_heatmap, processes, file_extension = "png") {
  # Construct file name
  combined_name <- paste(processes, collapse = "_")
  file_path <- file.path(output_dir_heatmap, paste("heatmap_", combined_name, ".", file_extension, sep = ""))
  
  # Print file path to verify it is being created correctly
  print(paste("Saving file to:", file_path))
  
  # Open a graphics device, e.g., png
  if(file_extension == "png") {
    png(filename = file_path, width = 10, height = 15, units = "in", res = 300, bg = "transparent")
  } else if (file_extension == "pdf") {
    pdf(file_path, width = 10, height = 15, bg = "transparent")
  } else if (file_extension == "tiff") {
    tiff(filename = file_path, width = 10, height = 15, units = "in", res = 300, bg = "transparent")
  } else {
    stop("Unsupported file extension: ", file_extension)
  }
  
  # Render and save the heatmap
  draw(heatmap)
  
  # Close the device
  dev.off()
}
