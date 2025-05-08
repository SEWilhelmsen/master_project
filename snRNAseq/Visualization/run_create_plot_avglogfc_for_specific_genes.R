# Create plot of average log fold change for multiple genes based on combined data 

# Source necessary scripts
#########################################################################
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/load_libraries_run_create_barplot_go_process.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/define_go_terms_and_find_genes.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/find_specific_genes.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/create_plot_for_specific_genes_for_combined_plot.R")
# source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_t_test_for_go_process.R")
# source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_anova_tukey_for_go_process.R")
# source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/create_barplot_go_process.R")
# source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/process_seurat_data.R")

# Set parameters and prepare data
############################################################################
go_process_of_interest <- 'tca'  # Change process 

output_dir_plot <- "C:/Users/siljeew/Master_project/snRNAseq/Plots/Gene_process"
file_path_to_data_avg_logfc <- "C:/Users/siljeew/Master_project/snRNAseq/Data/mouse_vcm_all_genes_avg_logFC.csv"
data_avg_logfc <- read_csv2(file_path_to_data_avg_logfc)
# Rename the first column using its position
colnames(data_avg_logfc)[1] <- "gene"
head(data_avg_logfc)

# Subset the data to include only the process genes
go_process_of_interest_genes <- gene_list[[go_process_of_interest]]

# Validate if the process of interest exist in the gene_list
if (is.null(go_process_of_interest_genes)) {
  stop(paste("GO term not found for the process:", go_process_of_interest))
}
print(go_process_of_interest_genes)


# Validate if the genes are present in the dataset
genes_in_data <- go_process_of_interest_genes %in% data_avg_logfc$gene
if (!all(genes_in_data)) {
  missing_genes <- go_process_of_interest_genes[!genes_in_data]
  stop(paste("The following genes are missing in the data:", paste(missing_genes, collapse = ", ")))
}

View(data_avg_logfc)
colnames(data_avg_logfc)

# Rename columns
data_avg_logfc <- data_avg_logfc %>%
  dplyr::rename(
    "6 Hours" = logFC6_hours,
    "12 Hours" = logFC12_hours,
    "1 Day" = logFC1_day,
    "3 Days" = logFC3_days,
    "1 Week" = logFC1_week,
    "3 Weeks" = logFC3_weeks
  )

View(data_avg_logfc)


# Create plot
#################################################################
# # Open the script for customisation: 
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/create_plot_for_genes_avglogFC.R")


# go_genes_plot <- create_go_genes_plot(data_avg_logfc, go_process_of_interest_genes)
# print(go_genes_plot)
# ggsave(file.path(output_dir_plot, paste(go_process_of_interest, "_genes_expression.png", sep = "")), plot = go_genes_plot, width = 8, height = 6)
# 
# ggsave(file.path(output_dir_plot, paste(go_process_of_interest, "_genes_expression.pdf", sep = "")), plot = go_genes_plot, width = 8, height = 6)


### Combine the two plots and save
#########################################################################
# combined_plot <- plot_grid(go_genes_plot, go_plot, align = "v", ncol = 1, rel_heights = c(1,1))
# combined_plot
# 
# 
# ggsave(file.path(output_dir_plot, paste(go_process_of_interest, "_and_genes_over_time.png", sep = "")), plot = combined_plot, width = 10, height = 14)
# ggsave(file.path(output_dir_plot, paste(go_process_of_interest, "_and_genes_over_time.pdf", sep = "")), plot = combined_plot, width = 10, height = 14)
