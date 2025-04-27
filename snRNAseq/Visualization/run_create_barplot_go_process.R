# Create a bar plot for one go process
# Silje Wilhelmsen


# Source necessary scripts
#########################################################################
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/load_libraries_run_create_barplot_go_process.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/create_barplot_go_process.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/process_seurat_data.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/define_go_terms_and_find_genes.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_t_test_for_go_process.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_anova_tukey_for_go_process.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/find_specific_genes.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/create_plot_for_specific_genes_for_combined_plot.R")

file.edit("C:/Users/siljeew/Master_project/snRNAseq/Analysis/define_go_terms_and_find_genes.R")
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/create_plot_for_specific_genes_for_combined_plot.R")



### Set parameters 
#########################################################################
go_process_of_interest <- 'fatty_acid_transport' # Change process of interest 
output_dir_plot <- "C:/Users/siljeew/Master_project/snRNAseq/Plots/Gene_process"

# Ensure output plot directory exists
if (!dir.exists(output_dir_plot)) {
  dir.create(output_dir_plot, recursive = TRUE)
}

# Define file paths for Seurat objects
file_paths <- list(
  "C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_6h_vcm_with_stress_status.Rds",
  "C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_12h_vcm_with_stress_status.Rds",
  "C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_1d_vcm_with_stress_status.Rds",
  "C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_3d_vcm_with_stress_status.Rds",
  "C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_1w_vcm_with_stress_status.Rds",
  "C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_3w_vcm_with_stress_status.Rds"
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


# Process and combine data from all Seurat objects
combined_data_long <- do.call(rbind, lapply(file_paths, process_seurat_object, genes_of_interest = genes_of_interest))



### Perform t-test to retrieve means
########################################################################
# Ensure the output directory exists for the data
output_dir_data <- "C:/Users/siljeew/Master_project/snRNAseq/Data"
if (!dir.exists(output_dir_data)) {
  dir.create(output_dir_data, recursive = TRUE)
}

# Set parameters
output_csv_path <- file.path(output_dir_data, "t_test_go_process_results.csv")
output_xlsx_path <- file.path(output_dir_data, "t_test_go_process_results.xlsx")

# Perform the t-tests with Benjamini-Hochberg adjustment
t_test_results <- perform_t_test(
  data = combined_data_long,
  group_col = "condition",
  value_col = "expression",
  time_col = "time_point",
  go_process_of_interest = go_process_of_interest,
  output_csv_path = output_csv_path,
  output_xlsx_path = output_xlsx_path
)

# Inspect results
t_test_go_process <- t_test_results %>%
  filter(process == go_process_of_interest)

print(t_test_go_process)

# Add a new column for significance label
t_test_go_process <- t_test_go_process %>%
  mutate(significance_label = ifelse(`p_adj` < 0.001, "***", 
                                     ifelse(`p_adj` < 0.01, "**", 
                                            ifelse(`p_adj` < 0.05, "*", "ns"))))
View(t_test_go_process)



### Perform ANOVA + TukeyHSD for mean comparison
########################################################################
# Define paths for saving results
anova_csv_path <- file.path(output_dir_data, "anova_go_process_results.csv")
anova_xlsx_path <- file.path(output_dir_data, "anova_go_process_results.xlsx")
tukey_csv_path <- file.path(output_dir_data, "tukey_go_process_results.csv")
tukey_xlsx_path <- file.path(output_dir_data, "tukey_go_process_results.xlsx")


# Perform the test
anova_tukey_process_result <- perform_anova_tukey_for_go_process(
  data = combined_data_long,  # Ensure this data is formatted correctly
  group_col = "condition",
  value_col = "expression",
  time_col = "time_point",
  go_process_of_interest = go_process_of_interest,
  anova_csv_path = anova_csv_path, 
  anova_xlsx_path = anova_xlsx_path, 
  tukey_csv_path = tukey_csv_path, 
  tukey_xlsx_path = tukey_xlsx_path
)




### Process Seurat objects
#########################################################################
# Aggregate gene expressions by time point and condition
go_process_aggregated <- combined_data_long %>%
  group_by(gene, time_point, condition) %>%
  summarize(mean_expression = mean(expression, na.rm = TRUE)) %>%
  ungroup()

# Summarize overall level by averaging gene expression levels
go_process_summary <- go_process_aggregated %>%
  group_by(time_point, condition) %>%
  summarize(mean_expression = mean(mean_expression, na.rm = TRUE), .groups = 'drop') %>%
  mutate(process = go_process_of_interest)

# Ensure 'condition' and time_point is correctly set as a factor
go_process_summary$time_point <- factor(go_process_summary$time_point, levels = unique(go_process_summary$time_point))
go_process_summary$condition <- factor(go_process_summary$condition, levels = c("SHAM", "ORAB"))

View(go_process_summary)


# Use p.adj from t-test
####################################################################
View(t_test_go_process)

p_adj_value <- t_test_go_process %>%
  filter(process == go_process_of_interest) %>%
  dplyr::select(process, time_point, p_adj, significance_label)

go_process_summary <- go_process_summary %>%
  left_join(p_adj_value, by = c("time_point", "process"))

View(go_process_summary)

# Go to create plot:
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/create_barplot_go_process.R")




# # Combine ANOVA and Tukey results with mean expression
# ###############################################################################
# View(go_process_summary)
# 
# # Convert factors to characters 
# anova_tukey_process_result <- anova_tukey_process_result %>%
#   mutate(Timepoint = as.character(Timepoint),
#          group1 = as.character(group1),
#          group2 = as.character(group2))
# 
# go_process_summary <- go_process_summary %>%
#   mutate(time_point = as.character(time_point))
# 
# view(anova_tukey_process_result)
# 
# # Insert mean expression from go_process_summary
# source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/insert_mean_in_tukey_result_go_process.R")
# # file.edit("C:/Users/siljeew/Master_project/snRNAseq/Analysis/insert_mean_in_tukey_result_go_process.R")
# 
# # Add a new column for significance label
# anova_tukey_process_result <- anova_tukey_process_result %>%
#   mutate(significance_label = ifelse(`p adj` < 0.001, "***", 
#                                             ifelse(`p adj` < 0.01, "**", 
#                                                    ifelse(`p adj` < 0.05, "*", "ns"))))
# 
# # Inspect columns
# view(anova_tukey_process_result)
# 
# # Insert p adj and significance label from anova_tukey_process_result
# source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/insert_significance_label_to_go_process_summary.R")
# 
# file.edit("C:/Users/siljeew/Master_project/snRNAseq/Analysis/insert_significance_label_to_go_process_summary.R")


### Generate and save process plot
#########################################################################
# Open the script for customisation: 
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/create_barplot_go_process.R")


# go_plot <- generate_plot(go_process_summary, p_values_for_plot, go_process_of_interest)
# print(go_plot)
# ggsave(file.path(output_dir_plot, paste(go_process_of_interest, "_levels_over_time_by_group.png", sep = "")), plot = go_plot, width = 8, height = 6)
# ggsave(file.path(output_dir_plot, paste(go_process_of_interest, "_levels_over_time_by_group.pdf", sep = "")), plot = go_plot, width = 8, height = 6)



