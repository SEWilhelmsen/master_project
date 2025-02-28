# Create bar plot of one gene ontology process
# Silje Wilhelmsen

library(dplyr)
library(tidyr)
library(Seurat)
library(ggpubr)
library(ggplot2)
library(tibble)


# # Function to adjust Seurat metadata
# process_metadata <- function(seurat_object) {
#   seurat_object@meta.data <- seurat_object@meta.data %>%
#     mutate(time_point = sub(" - .*", "", orig.ident), 
#            condition = sub(".*- ", "", orig.ident))
#   return(seurat_object)
# }
# 
# 
# # Function to extract and filter gene expression data
# process_expression_data <- function(seurat_object, genes_of_interest) {
#   data <- GetAssayData(seurat_object, slot = "data") %>%
#     as.data.frame()
#   
#   # Keep only the genes from genes_of_interest
#   data <- data[rownames(data) %in% toupper(genes_of_interest), ]
#   print(paste("Number of genes in filtered data:", nrow(data)))
#   
#   rownames(data) <- toupper(rownames(data))
#   
#   # Convert data to long format
#   data_long <- data %>%
#     rownames_to_column(var = "gene") %>%
#     pivot_longer(cols = -gene, names_to = "cell", values_to = "expression") %>%
#     mutate(time_point = seurat_object@meta.data$time_point[match(cell, rownames(seurat_object@meta.data))],
#            condition = seurat_object@meta.data$condition[match(cell, rownames(seurat_object@meta.data))])
#   
#   return(data_long)
# }
# 
# # Process each Seurat object and return combined data
# process_seurat_object <- function(file_path, genes_of_interest) {
#   seurat_object <- readRDS(file_path)
#   seurat_object <- process_metadata(seurat_object)
#   data_long <- process_expression_data(seurat_object, genes_of_interest)
#   return(data_long)
# }
# 
# # Combine data from all Seurat objects
# combined_data_long <- do.call(rbind, lapply(file_paths, process_seurat_object, genes_of_interest = genes_of_interest))


# Function to generate plots with significance stars
generate_plot <- function(summary_data, significance_data, go_process_of_interest) {
  significance_data <- significance_data %>%
    filter(process == go_process_of_interest) %>%
    dplyr::select(time_point, group1, group2, p_value) %>%
    mutate(
      significance = case_when(
        p_value < 0.001 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  significance_data <- significance_data %>%
    left_join(
      summary_data %>%
        group_by(time_point) %>%
        summarize(max_mean_expression = max(mean_expression, na.rm = TRUE)),
      by = "time_point"
    ) %>%
    mutate(
      y.position = max_mean_expression * 1.05
    ) %>%
    dplyr::select(time_point, group1, group2, p_value, significance, y.position)
  
  summary_data$time_point <- factor(summary_data$time_point,
                                    levels = c("6 Hours",
                                               "12 Hours",
                                               "1 Day",
                                               "3 Days",
                                               "1 week",
                                               "3 Weeks"))
  
  max_y_position <- max(significance_data$y.position, na.rm = TRUE)
  
  go_plot <- ggplot(summary_data, aes(x = time_point, y = mean_expression, fill = condition)) +
    geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
    labs(x = NULL, y = paste(go_process_of_interest, "Level")) +
    theme_minimal() +
    scale_fill_manual(values = c("SHAM" = "grey", "ORAB" = "coral")) +
    scale_y_continuous(limits = c(0, max_y_position * 1.1)) +  # Add extra space above the bars
    geom_text(
      data = significance_data,
      aes(x = time_point, y = y.position, label = significance, fill = NULL),  # Remove fill from text layer
      vjust = -0.5,
      position = position_dodge(width = 0.75)) +
    scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
  
  return(go_plot)
}
