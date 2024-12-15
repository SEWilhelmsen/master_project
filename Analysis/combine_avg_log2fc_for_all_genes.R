## Combine all genes from .xlsx, only the column:  avg_log2FC (*ikke ferdig*)
# Silje Wilhelmsen

# Load necessary libraries
library(readxl)
library(dplyr)
library(tibble)

# Function to read data from a list of filenames and time points
read_mouse_data <- function(file_paths, time_points) {
  data_list <- list()
  
  for (i in seq_along(file_paths)) {
    if (file.exists(file_paths[i])) {
      data_list[[i]] <- readxl::read_excel(file_paths[i])
    } else {
      warning(paste("File does not exist:", file_paths[i]))
    }
  }
  
  # Ensure the data_list lengths match
  if (length(data_list) == length(time_points)) {
    names(data_list) <- time_points
  } else {
    stop("Lengths of data_list and time_points do not match")
  }
  
  return(data_list)
}



# Function to normalize gene names to uppercase and rename columns
normalize_and_rename <- function(mouse_data_list) {
  renamed_data_list <- list()
  
  for (time_point in names(mouse_data_list)) {
    data <- mouse_data_list[[time_point]]
    
    # Normalize gene names to uppercase
    data <- data %>%
      mutate(gene = toupper(gene)) %>%
      dplyr::select(gene, avg_log2FC) %>%
      rename_at(vars(avg_log2FC), ~paste0("avg_log2fc_", time_point))
    
    renamed_data_list[[time_point]] <- data
  }
  
  return(renamed_data_list)
}




# # Import data
# #################
# 
# # mouse_6h_vcm_markers_all_genes <- read_excel("C:/Users/siljeew/snRNAseq/Data/mouse_6h_vcm_markers_all_genes.xlsx")
# # mouse_12h_vcm_markers_all_genes <- read_excel("C:/Users/siljeew/snRNAseq/Data/mouse_12h_vcm_markers_all_genes.xlsx")
# # mouse_1d_vcm_markers_all_genes <- read_excel("C:/Users/siljeew/snRNAseq/Data/mouse_1d_vcm_markers_all_genes.xlsx")
# mouse_3d_vcm_markers_all_genes <- read_excel("C:/Users/siljeew/snRNAseq/Data/mouse_3d_vcm_markers_all_genes.xlsx")
# mouse_1w_vcm_markers_all_genes <- read_excel("C:/Users/siljeew/snRNAseq/Data/mouse_1w_vcm_markers_all_genes.xlsx")
# mouse_3w_vcm_markers_all_genes <- read_excel("C:/Users/siljeew/snRNAseq/Data/mouse_3w_vcm_markers_all_genes.xlsx")
# 
# # Put the data into a list to use in the functions
# mouse_data_list <- list(
#   # '6_hours' = mouse_6h_vcm_markers_all_genes,
#   # '12_hours' = mouse_12h_vcm_markers_all_genes,
#   # '1_day' = mouse_1d_vcm_markers_all_genes,
#   '3_days' = mouse_3d_vcm_markers_all_genes,
#   '1_week' = mouse_1w_vcm_markers_all_genes,
#   '3_weeks' = mouse_3w_vcm_markers_all_genes
# )
# 
# 
# 
# # Normalize gene names to uppercase
# #######################################
# normalize_and_rename <- function(mouse_data_list){
#   # Initialize a list to store the renamed data frames
#   renamed_data_list <- list()
#   
#   for (time_point in names(mouse_data_list)){
#     data <- mouse_data_list[[time_point]]
#     
#     # Normalize gene names to uppercase
#     data$gene <- toupper(data$gene)
#     
#     # Select and rename columns
#     data_avg_log2fc <- data[, c("gene", "avg_log2FC")]
#     colnames(data_avg_log2fc)[2] <- paste0("avg_log2fc_", time_point)
#     
#     # Add the renames data frame to the list
#     renamed_data_list[[time_point]] <- data_avg_log2fc
#     
#   }
#   
#   # Return the list of renamed data frames 
#   return(renamed_data_list)
# }
# 
# normalized_data <- normalize_and_rename(mouse_data_list)
#   
# merged_data <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), normalized_data)  
# View(merged_data)
# 
# # Convert the "gene" column to row names
# merged_data <- merged_data %>%
#   column_to_rownames(var = "gene")
# 
# # # Normalize gene names to uppercase
# # mouse_3d_vcm_markers_all_genes$gene <- toupper(mouse_3d_vcm_markers_all_genes$gene)
# # mouse_1w_vcm_markers_all_genes$gene <- toupper(mouse_1w_vcm_markers_all_genes$gene)
# # mouse_3w_vcm_markers_all_genes$gene <- toupper(mouse_3w_vcm_markers_all_genes$gene)
# # 
# # # Use base R to select and rename columns
# # # For 3 days
# # mouse_3d_vcm_all_genes_avg_log2fc <- mouse_3d_vcm_markers_all_genes[, c("gene", "avg_log2FC")]
# # colnames(mouse_3d_vcm_all_genes_avg_log2fc)[2] <- "avg_log2fc_3d"
# # 
# # # For 1 week
# # mouse_1w_vcm_all_genes_avg_log2fc <- mouse_1w_vcm_markers_all_genes[, c("gene", "avg_log2FC")]
# # colnames(mouse_1w_vcm_all_genes_avg_log2fc)[2] <- "avg_log2fc_1w"
# # 
# # # For 3 weeks
# # mouse_3w_vcm_all_genes_avg_log2fc <- mouse_3w_vcm_markers_all_genes[, c("gene", "avg_log2FC")]
# # colnames(mouse_3w_vcm_all_genes_avg_log2fc)[2] <- "avg_log2fc_3w"
# 
# 
# # # Check if columns are correctly renamed
# # print(head(mouse_3d_vcm_all_genes_avg_log2fc))
# # print(head(mouse_1w_vcm_all_genes_avg_log2fc))
# # print(head(mouse_3w_vcm_all_genes_avg_log2fc))
# 
# 
# # Merge the data based on gene names
# # mouse_vcm_all_genes_avg_log2fc <- merge(
# #   merge(mouse_3d_vcm_all_genes_avg_log2fc, mouse_1w_vcm_all_genes_avg_log2fc, by = "gene", all = TRUE),
# #   mouse_3w_vcm_all_genes_avg_log2fc,
# #   by = "gene",
# #   all = TRUE
# # )
# 
# # # Print the merged data to verify the merge
# # print(head(mouse_vcm_all_genes_avg_log2fc))
# 
# # # Convert the "gene" column to row names
# # mouse_vcm_all_genes_avg_log2fc <- mouse_vcm_all_genes_avg_log2fc %>%
# #   column_to_rownames(var = "gene")
# 
# # # Save the merged data to a CSV file
# # write.csv(mouse_vcm_all_genes_avg_log2fc, "mouse_vcm_all_genes_avg_log2fc.csv", row.names = FALSE)
# # 
# # # Save the merged data to an RDS file for efficient storage and retrieval in R
# # saveRDS(mouse_vcm_all_genes_avg_log2fc, "mouse_vcm_all_genes_avg_log2fc.rds")
# 
