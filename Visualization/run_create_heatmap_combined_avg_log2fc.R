# Run create heatmap 
# Silje Wilhelmsen

# Source functions
########################################################################
source("C:/Users/siljeew/snRNAseq/Visualization/load_libraries_run_create_heatmap_combined_avg_log2fc.R")
source("C:/Users/siljeew/snRNAseq/Analysis/define_go_terms_and_find_genes.R")
source("C:/Users/siljeew/snRNAseq/Visualization/gene_process_mapping.R")
source("C:/Users/siljeew/snRNAseq/Visualization/create_heatmap_combined_avg_log2fc.R")

### Set parameters
#######################################################################
output_dir_plot <- "C:/Users/siljeew/snRNAseq/Plots"

# Ensure the directory exists
if (!dir.exists(output_dir_plot)) {
  dir.create(output_dir_plot, recursive = TRUE)
}

#Load data
mouse_vcm_all_genes_avg_log2fc <- readRDS("C:/Users/siljeew/snRNAseq/tmp/mouse_vcm_all_genes_avg_log2fc.rds")
head(mouse_vcm_all_genes_avg_log2fc)

# Define go processes of interest
go_processes_of_interest_for_heatmap <- c('glycolysis', 'tca')

# Retrieve genes of interest
genes_of_interest_for_heatmap_vector <- get_genes_of_interest(go_processes_of_interest_for_heatmap, define_go_process, get_genes_for_go_process)
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
  gene <- toupper(gene)  # Ensure consistent uppercase case
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

# Define a sorting order
sorting_order <- factor(go_association_vector, levels = c("glycolysis", "glycolysis,tca", "tca"))

# Order the matrix based on the sorting vector
ordered_indices <- order(sorting_order)

# Reorder heatmap matrix and process vector
heatmap_matrix_numeric <- heatmap_matrix_numeric[ordered_indices, ]
go_association_vector <- go_association_vector[ordered_indices]

# Verify sorting
print(head(rownames(heatmap_matrix_numeric)))
print(head(go_association_vector))

### Create and save heatmap
#######################################################################
heatmap_plot <- create_heatmap(heatmap_matrix_numeric, go_association_vector)

ggsave(file.path(output_dir_plot, paste(go_processes_of_interest_for_heatmap, "_heatmap.png", sep = "")), plot = combined_plot, width = 10, height = 7)

