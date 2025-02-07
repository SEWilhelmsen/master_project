# Run create heatmap for transcription factors 
# Silje Wilhelmsen

# Source functions
########################################################################
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/load_libraries_run_create_heatmap_combined_avg_log2fc.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/define_go_terms_and_find_genes.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/gene_process_mapping.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/create_heatmap_combined_avg_log2fc_transcription_factors.R")

### Set parameters
#######################################################################
output_dir_heatmap <- "C:/Users/siljeew/Master_project/snRNAseq/Plots/Heatmap"

# Ensure the directory exists
if (!dir.exists(output_dir_heatmap)) {
  dir.create(output_dir_heatmap, recursive = TRUE)
}

#Load data
mouse_vcm_all_genes_avg_log2fc <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_vcm_all_genes_avg_log2fc.rds")
head(mouse_vcm_all_genes_avg_log2fc)

# Define go processes of interest
go_process_of_interest_for_heatmap <- 'transcription_coregulator_activity'


# Retrieve genes of interest
genes_of_interest_for_heatmap_vector <- get_genes_of_interest(go_process_of_interest_for_heatmap, define_go_process, get_genes_for_go_process)
View(genes_of_interest_for_heatmap_vector)

# Process `genes_of_interest_for_heatmap_vector` to ensure it's flat
genes_of_interest_for_heatmap_vector_unlisted <- unlist(genes_of_interest_for_heatmap_vector, use.names = FALSE)

# Create gene process mapping using the retrieved genes
gene_process_mapping <- create_gene_process_mapping(genes_of_interest_for_heatmap_vector)


# Prepare the heatmap matrix
heatmap_matrix_numeric <- prepare_heatmap_matrix(mouse_vcm_all_genes_avg_log2fc, genes_of_interest_for_heatmap_vector_unlisted)
nrow(heatmap_matrix_numeric)

# Check the column names and change to shorter and more readable names
colnames(heatmap_matrix_numeric) <- c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")




# Create a vector for gene ontology mapping
#######################################################################
# Must be executed post `gene_process_mapping` calculation
go_association_vector <- sapply(rownames(heatmap_matrix_numeric), function(gene) {
  gene <- toupper(gene)  
  if (gene %in% names(gene_process_mapping)) {
    processes <- gene_process_mapping[[gene]]
    unique_processes <- unique(processes)
    clean_process_name <- paste(unique_processes, collapse = ",")
    clean_process_name <- gsub("[0-9]+", "", clean_process_name)
    return(clean_process_name)
  } else {
    print(paste("No process found for gene:", gene))
    return("None")
  }
}, USE.NAMES = FALSE)




# # Verify vector
# print(head(rownames(heatmap_matrix_numeric)))
# print(head(go_association_vector))




# Filter data (optional)
###################################################################33
# Filter rows with less than four NA values
rows_to_keep <- rowSums(is.na(heatmap_matrix_numeric)) < 4
filtered_heatmap_matrix <- heatmap_matrix_numeric[rows_to_keep, ]
go_association_vector_filtered <- go_association_vector[rows_to_keep]

# Exclude rows with NA in the first and last columns
rows_to_exclude <- is.na(filtered_heatmap_matrix[, 1]) & is.na(filtered_heatmap_matrix[, ncol(filtered_heatmap_matrix)])
filtered_matrix <- filtered_heatmap_matrix[!rows_to_exclude, ]
go_association_vector_filtered <- go_association_vector_filtered[!rows_to_exclude]

# Print dimensions to check results
cat("Original dimensions:", dim(heatmap_matrix_numeric), "\n")
cat("After filtering for NAs in fewer than four columns:", dim(filtered_heatmap_matrix), "\n")
cat("Filtered dimensions after excluding rows with NA in first and last column:", dim(filtered_matrix), "\n")

# Check lengths match
stopifnot(nrow(filtered_matrix) == length(go_association_vector_filtered))





### Create and save heatmap (filtered)
#######################################################################
# Create heatmap
heatmap_plot <- create_heatmap(filtered_matrix, go_association_vector_filtered)
print(heatmap_plot)


# Save heatmap
save_heatmap_complex(filtered_matrix, heatmap_plot, output_dir_heatmap, go_process_of_interest_for_heatmap, file_extension = "png")



# ### Create and save heatmap (not filtered)
# #######################################################################
# # Create heatmap
# heatmap_plot <- create_heatmap(heatmap_matrix_numeric, go_association_vector)
# print(heatmap_plot)
# 
# 
# # Save heatmap
# save_heatmap_complex(heatmap_matrix_numeric, heatmap_plot, output_dir_heatmap, go_process_of_interest_for_heatmap, file_extension = "png")
# 
