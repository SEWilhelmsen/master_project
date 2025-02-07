## Create heatmap from combined_avg_log2fc for transcription factors
# Silje Wilhelmsen



# Prepare heatmap matrix should consider global variable dependency
prepare_heatmap_matrix <- function(data, genes_of_interest_for_heatmap_vector) {
  # Ensure both lists are uppercase
  genes_of_interest_for_heatmap_vector <- toupper(genes_of_interest_for_heatmap_vector)

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
create_heatmap <- function(filtered_matrix, go_association_vector_filtered) {
  if (!is.matrix(filtered_matrix)) {
    stop("Input data must be a matrix.")
  }
  
  if (nrow(filtered_matrix) != length(go_association_vector_filtered)) {
    stop("Mismatch between number of matrix rows and association vector length.")
  }
  
  # Define color function
  col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  # Determine height based on number of rows
  num_rows <- nrow(filtered_matrix)
  plot_height <- unit(num_rows * 5, "mm") # Adjust multiplier as needed
  
  
  # Create the heatmap
  heatmap <- Heatmap(
    filtered_matrix,
    name = "log2 fold change",
    col = col_fun,
    na_col = "grey",
    cluster_rows = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 10),
    row_gap = unit(10, "mm"),
    cluster_columns = FALSE,
    show_column_names = TRUE,
    column_names_centered = TRUE,
    column_names_rot = 30,
    column_title = "Differential Gene Expression ORAB vs SHAM",
    show_column_dend = FALSE,
    width = unit(15, "cm"),
    height = plot_height,
    heatmap_legend_param = list(
      title = "average log 2 fold change",
      legend_direction = "vertical",
      legend_height = unit(4, "cm"),
      at = c(-1, -0.5, 0, 0.5, 1),
      labels = c("-1", "Higher expression in SHAM",  "0", "Higher expression in ORAB", "1")
    ),
    gap = unit(1, "mm")
  )
  
  # Create custom legend for NA
  na_legend <- Legend(
    labels = "NA",
    title = "",
    legend_gp = gpar(fill = "grey"),
    direction = "vertical"
  )
  
  
  # Draw heatmap with legends
  draw(heatmap, annotation_legend_list = list(na_legend), merge_legend = TRUE)
}




# Function to save a heatmap
save_heatmap_complex <- function(filtered_matrix, heatmap, output_dir_heatmap, go_process_of_interest_for_heatmap, file_extension = "png") {
  
  # Calculate height based on number of rows
  num_rows <- nrow(filtered_matrix)
  plot_height <- num_rows * (5 / 25.4) # Convert mm to inches for the device
  
  # Construct file name
  combined_name <- paste(go_process_of_interest_for_heatmap, num_rows, collapse = "_", sep = "")
  file_path <- file.path(output_dir_heatmap, paste("heatmap_", combined_name, ".", file_extension, sep = ""))
  
  # Print file path to verify it is being created correctly
  print(paste("Saving file to:", file_path))

  
  # Open a graphics device, e.g., png
  if(file_extension == "png") {
    png(filename = file_path, width = 15, height = plot_height, units = "in", res = 300, bg = "transparent")
  } else if (file_extension == "pdf") {
    pdf(file_path, width = 15, height = plot_height, bg = "transparent")
  } else if (file_extension == "tiff") {
    tiff(filename = file_path, width = 15, height = plot_height, units = "", res = 300, bg = "transparent")
  } else {
    stop("Unsupported file extension: ", file_extension)
  }
  
  # Render and save the heatmap
  draw(heatmap)
  
  # Close the device
  dev.off()
}

