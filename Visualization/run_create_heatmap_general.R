# Run create heatmap general
# Silje Wilhelmsen

source("C:/Users/siljeew/snRNAseq/Visualization/load_libraries_run_create_heatmap_general.R")
source("C:/Users/siljeew/snRNAseq/Analysis/combine_avg_log2fc_for_all_genes.R")
source("C:/Users/siljeew/snRNAseq/Visualization/create_heatmap_general.R")

### Define GO terms and find related genes 
source("C:/Users/siljeew/snRNAseq/Analysis/define_go_terms_and_find_genes.R")

# For finding specific genes
source("C:/Users/siljeew/snRNAseq/Analysis/find_specific_genes.R")

### Set parameters
#######################################################################
output_dir_plot <- "C:/Users/siljeew/snRNAseq/Plots"

# Ensure the directory exists
if (!dir.exists(output_dir_plot)) {
  dir.create(output_dir_plot, recursive = TRUE)
}


#### Prepare data 
#######################################################################

# File paths and corresponding time points
file_paths <- c(
  "C:/Users/siljeew/snRNAseq/Data/mouse_3d_vcm_markers_all_genes.xlsx",
  "C:/Users/siljeew/snRNAseq/Data/mouse_1w_vcm_markers_all_genes.xlsx",
  "C:/Users/siljeew/snRNAseq/Data/mouse_3w_vcm_markers_all_genes.xlsx"
)

time_points <- c("3_days", "1_week", "3_weeks")

# Read the mouse data
mouse_data_list <- read_mouse_data(file_paths, time_points)


# Normalize and rename the data
normalized_data <- normalize_and_rename(mouse_data_list)

# Merge all data frames in the list into one based on gene names
merged_data <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), normalized_data)

# View the merged data
View(merged_data)

# Convert the "gene" column to row names
merged_data %>% column_to_rownames(var = "gene")
print(head(merged_data))


### Find genes for multiple gene ontology processes
#######################################################################
go_terms_of_interest_list <- c('tca', 'glycolysis', 'pdh', 'etc') # How to choose multiple terms in the list?

go_genes_df <- data.frame(gene = character(), process = character(), stringsAsFactors = FALSE)

if (is.null(go_term)) {
  stop(paste("GO term not found for the process:", go_process_of_interest))
}

# Loop through each GO process, retrieve genes, and add them to the data frame
for (term in go_terms_list_of_interest) {
  go_term <- define_go_process[[term]]
  if (is.null(go_term)) {
    stop(paste("GO term not found for the process:", term))
  }
  genes <- get_genes_for_go_process(go_term)
  process_genes_df <- data.frame(gene = genes, process = term, stringsAsFactors = FALSE)
  go_genes_df <- rbind(go_genes_df, process_genes_df)
}

# Print the data frame
head(go_genes_df)


### Filter the data to only include genes in the gene ontology processes
#########################################################################

# Filter the genes in go_genes_df
filtered_data <- filter_and_assign_process(merged_data, go_genes_df)
View(filtered_data)

### OR: Filter the data to only include specific genes
########################################################################

regulated_genes_df <- data.frame(gene = character(), process = character(), stringsAsFactors = FALSE)

if (is.null(go_term)) {
  stop(paste("GO term not found for the process:", go_process_of_interest))
}

# Loop through each GO process, retrieve genes, and add them to the data frame
for (term in go_terms_list_of_interest) {
  go_term <- define_go_process[[term]]
  if (is.null(go_term)) {
    stop(paste("GO term not found for the process:", term))
  }
  genes <- get_genes_for_go_process(go_term)
  process_genes_df <- data.frame(gene = genes, process = term, stringsAsFactors = FALSE)
  go_genes_df <- rbind(go_genes_df, process_genes_df)
}

# Print the data frame
head(go_genes_df)



# Subset the data to include only the process genes
go_process_of_interest_genes <- gene_list[[go_terms_of_interest_list]]

# Validate if the process of interest exist in the gene_list
if (is.null(go_process_of_interest_genes)) {
  stop(paste("GO term not found for the process:", go_process_of_interest))
}
print(go_process_of_interest_genes)

# Validate if the genes are present in the dataset
# data <- read.csv(file_path_to_data_avg_log2fc)
# genes_of_interest <- go_process_of_interest_genes


genes_in_data <- go_process_of_interest_genes %in% data$gene
if (!all(genes_in_data)) {
  missing_genes <- go_process_of_interest_genes[!genes_in_data]
  stop(paste("The following genes are missing in the data:", paste(missing_genes, collapse = ", ")))
}



### Create and save heatmap
#######################################################################
# Prepare the gene matrix for the heatmap
gene_matrix_ordered <- prepare_gene_matrix(filtered_data)
print(colnames(gene_matrix_ordered)) # Confirm the column order
print(head(gene_matrix_ordered))

# Create and print the heatmap
heatmap_plot <- create_heatmap(gene_matrix_ordered, filtered_data$go_category)
print(heatmap_plot)


ggsave(file.path(output_dir_plot, paste("heatmap.png", sep = "")), plot = plot1, width = 8, height = 6)
