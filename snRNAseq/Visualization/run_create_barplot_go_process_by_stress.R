# Create a bar plot for one go process
# Silje Wilhelmsen

# Somewhere the names are capital and small letters, wtf? 

# Source necessary scripts
#########################################################################
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/load_libraries_run_create_barplot_go_process.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/create_barplot_go_process.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/process_seurat_data.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/define_go_terms_and_find_genes.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_t_test_for_go_process.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/find_specific_genes.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/create_plot_for_specific_genes_for_combined_plot.R")




### Set parameters 
#########################################################################
go_process_of_interest <- 'hydroxybutyrate_dehydrogenase_activity'
output_dir_plot <- "C:/Users/siljeew/Master_project/snRNAseq/Plots/Gene_process"

# Ensure output plot directory exists
if (!dir.exists(output_dir_plot)) {
  dir.create(output_dir_plot, recursive = TRUE)
}

# Define file paths for Seurat objects
file_paths <- list(
  "C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_6h_vcm.Rds",
  "C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_12h_vcm.Rds",
  "C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_1d_vcm.Rds",
  "C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_vcm_3d.Rds",
  "C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_vcm_1w.Rds",
  "C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_vcm_3w.Rds"
)


### Retrieve genes for the gene ontology process
#########################################################################
# Retrieve genes for the GO term of interest
go_term <- define_go_process[[go_process_of_interest]]
if (is.null(go_term)) {
  stop(paste("GO term not found for the process:", go_process_of_interest))
}

# Get genes related to the selected GO process
genes_of_interest <- get_genes_for_go_process(go_term)
print(paste("Genes of interest for process", go_process_of_interest, "are:", paste(genes_of_interest, collapse = ", ")))


#genes_of_interest <- c("MYH7", "NPPA", "NPPB", "ANKRD1", "HK1", "HK2", "GCK", "PFK1", "PFKM", "PK")

### Process seurat objects
#########################################################################
# Process and combine data from all Seurat objects
combined_data_long <- do.call(rbind, lapply(file_paths, process_seurat_object, genes_of_interest = genes_of_interest))



### Perform t-test
########################################################################
# Ensure the output directory exists for the data 
output_dir_data <- "C:/Users/siljeew/Master_project/snRNAseq/Data"
if (!dir.exists(output_dir_data)) {
  dir.create(output_dir_data, recursive = TRUE)
}

# Set parameters
output_csv_path <- file.path(output_dir_data, "t_test_go_process_results.csv")
output_xlsx_path <- file.path(output_dir_data, "t_test_go_process_results.xlsx")

# Perform the t-tests
p_values_for_plot <- perform_t_test(
  data = combined_data_long, 
  group_col = "condition", 
  value_col = "expression", 
  time_col = "time_point", 
  go_process_of_interest = go_process_of_interest,
  output_csv_path = output_csv_path,
  output_xlsx_path = output_xlsx_path
)

tail(p_values_for_plot, 6)


### Process seurat objects
#########################################################################
# Aggregate gene expressions by time point and condition
go_process_aggregated <- combined_data_long %>%
  group_by(gene, time_point, condition) %>%
  summarize(mean_expression = mean(expression, na.rm = TRUE)) %>%
  ungroup()

# Summarize overall level by averaging gene expression levels
go_process_summary <- go_process_aggregated %>%
  group_by(time_point, condition) %>%
  summarize(mean_expression = mean(mean_expression, na.rm = TRUE), .groups = 'drop')

# Ensure 'condition' is correctly set as a factor
go_process_summary$time_point <- factor(go_process_summary$time_point, levels = unique(go_process_summary$time_point))
go_process_summary$condition <- factor(go_process_summary$condition, levels = c("SHAM", "ORAB"))

View(go_process_summary)
View(p_values_for_plot)



### Generate and save process plot
#########################################################################
go_plot <- generate_plot(go_process_summary, p_values_for_plot, go_process_of_interest)
print(go_plot)
ggsave(file.path(output_dir_plot, paste(go_process_of_interest, "_levels_over_time_by_group.png", sep = "")), plot = go_plot, width = 8, height = 6)
ggsave(file.path(output_dir_plot, paste(go_process_of_interest, "_levels_over_time_by_group.pdf", sep = "")), plot = go_plot, width = 8, height = 6)




### Create genes plot
#########################################################################
# Set parameters 
file_path_to_data_avg_log2fc <- "C:/Users/siljeew/Master_project/snRNAseq/Data/mouse_vcm_all_genes_avg_log2fc.csv"

data_avg_log2fc <- read_csv2(file_path_to_data_avg_log2fc)

# Rename the first column using its position
colnames(data_avg_log2fc)[1] <- "gene"
#View(data_avg_log2fc)


# Subset the data to include only the process genes
go_process_of_interest_genes <- gene_list[[go_process_of_interest]]

# Validate if the process of interest exist in the gene_list
if (is.null(go_process_of_interest_genes)) {
  stop(paste("GO term not found for the process:", go_process_of_interest))
}
print(go_process_of_interest_genes)


# Validate if the genes are present in the dataset
genes_in_data <- go_process_of_interest_genes %in% data_avg_log2fc$gene
if (!all(genes_in_data)) {
  missing_genes <- go_process_of_interest_genes[!genes_in_data]
  stop(paste("The following genes are missing in the data:", paste(missing_genes, collapse = ", ")))
}


# Create plot
go_genes_plot <- create_go_genes_plot(data_avg_log2fc, go_process_of_interest_genes)
print(go_genes_plot)
ggsave(file.path(output_dir_plot, paste(go_process_of_interest, "_genes_expression.png", sep = "")), plot = go_genes_plot, width = 8, height = 6)

ggsave(file.path(output_dir_plot, paste(go_process_of_interest, "_genes_expression.pdf", sep = "")), plot = go_genes_plot, width = 8, height = 6)


### Combine the two plots and save
#########################################################################
combined_plot <- plot_grid(go_genes_plot, go_plot, align = "v", ncol = 1, rel_heights = c(1,1))
combined_plot

ggsave(file.path(output_dir_plot, paste(go_process_of_interest, "_and_genes_over_time.png", sep = "")), plot = combined_plot, width = 8, height = 6)

ggsave(file.path(output_dir_plot, paste(go_process_of_interest, "_and_genes_over_time.pdf", sep = "")), plot = combined_plot, width = 8, height = 6)
