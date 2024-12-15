# Create heatmap of gene expression
# silje Wilhelmsen


# Load required libraries
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

# Function to filter and assign unique GO categories
filter_and_assign_process <- function(merged_data, go_genes_df) {
  # Remove duplicate genes by priority of process (in the order they are defined in define_go_process)
  go_genes_df <- go_genes_df %>% 
    dplyr::distinct(gene, .keep_all = TRUE)
  
  merged_data %>%
    filter(gene %in% go_genes_df$gene) %>%
    left_join(go_genes_df, by = "gene")
}

# Function to prepare gene matrix for the heatmap
prepare_gene_matrix <- function(data) {
  data_with_row_names <- data %>%
    column_to_rownames(var = "gene")
  
  as.matrix(data_with_row_names[, c("avg_log2fc_3_days", "avg_log2fc_1_week", "avg_log2fc_3_weeks")])
}

# Function to create heatmap
create_heatmap <- function(gene_matrix, go_categories) {
  Heatmap(gene_matrix, 
          name = "log2fc", 
          cluster_rows = FALSE, 
          show_row_names = TRUE,
          row_names_side = "left",
          row_split = go_categories,
          row_names_gp = gpar(fontsize = 6),
          row_gap = unit(2, "mm"),
          cluster_columns = FALSE,  
          show_column_names = TRUE,
          column_names_rot = 30,                        # Rotate column names 30 degrees
          column_title = "Expression AB vs SHAM",
          show_column_dend = FALSE,
          heatmap_legend_param = list(
            title = "average log2 fold change",
            legend_direction = "vertical",
            legend_height = unit(4, "cm")
          ),
          gap = unit(2, "mm"),
          heatmap_height = unit(0.3, "cm") * nrow(gene_matrix),
          col = colorRamp2(c(min(gene_matrix, na.rm = TRUE), 
                             0, 
                             max(gene_matrix, na.rm = TRUE)),
                           c("navy", "white", "red4")),
          na_col = "grey" # Color for missing values
  )
}

heatmap_plot <- create_heatmap(gene_matrix_ordered, filtered_data$go_category)
print(heatmap_plot)

# # Subset the data frame to include only the genes in all_genes_for_go_terms
# mouse_vcm_all_genes_avg_log2fc_filtered <- mouse_vcm_all_genes_avg_log2fc %>%
#   filter(gene %in% all_genes_for_go_terms)
# 
# # Create an identifier column for each GO category and determine their corresponding genes
# mouse_vcm_all_genes_avg_log2fc_filtered <- mouse_vcm_all_genes_avg_log2fc_filtered %>%
#   mutate(go_category = case_when(
#     gene %in% normalized_go_genes$glycolysis ~ 'glycolysis',
#     gene %in% normalized_go_genes$fao ~ 'fao',
#     gene %in% normalized_go_genes$tca ~ 'tca'
#   ))
# 
# print(head(mouse_vcm_all_genes_avg_log2fc_filtered)) # Check the first few rows of the filtered data frame
# 
# # Prepare the data for heatmap creation
# gene_matrix <- mouse_vcm_all_genes_avg_log2fc_filtered %>%
#   column_to_rownames(var = "gene")
# 
# # Convert the subset of data to a matrix
# gene_matrix_ordered <- as.matrix(gene_matrix[, c("avg_log2FC_24h", "avg_log2FC_3d", "avg_log2FC_1w", "avg_log2FC_3w")])
# 
# # Print the first few rows of the ordered gene matrix to confirm the column order
# print(colnames(gene_matrix_ordered))
# print(head(gene_matrix_ordered))
# 
# 
# 
# # Create a heatmap sorted by GO category
# create_heatmap <- Heatmap(gene_matrix_ordered, 
#                           name = "log2fc", 
#                           cluster_rows = FALSE, 
#                           show_row_names = TRUE,
#                           row_names_side = "left",
#                           row_split = mouse_vcm_all_genes_avg_log2fc_filtered$go_category,
#                           row_names_gp = gpar(fontsize = 12),
#                           row_gap = unit(2, "mm"),
#                           cluster_columns = FALSE,  
#                           show_column_names = TRUE,
#                           column_names_rot = 30,                          # Rotate column names 30 degrees
#                           column_title = "Expression AB vs SHAM",
#                           show_column_dend = FALSE,
#                           heatmap_legend_param = list(
#                             title = "avg_log2FC",
#                             legend_direction = "vertical",
#                             legend_height = unit(4, "cm")
#                           ),
#                           gap = unit(2, "mm"),
#                           heatmap_height = unit(5, "cm")*nrow(gene_matrix_ordered),
#                           col = colorRamp2(c(min(gene_matrix_ordered, na.rm = TRUE), 
#                                              0, 
#                                              max(gene_matrix_ordered, na.rm = TRUE)),
#                                            c("navy", "white", "red4")),
#                           na_col = "grey" # Color for missing values
# )
# 
# # Print the heatmap
# print(create_heatmap)
