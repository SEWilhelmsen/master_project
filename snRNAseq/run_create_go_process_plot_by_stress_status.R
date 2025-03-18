# Create a bar plot for go process grouped by stress status
# Silje Wilhelmsen

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
install.packages("ggpubr")
library(ggpubr)
library(openxlsx)

# Source necessary scripts
#########################################################################
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/load_libraries_run_create_barplot_go_process.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/define_go_terms_and_find_genes.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_anova_tukey_by_stress_status.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/find_specific_genes.R")



### Set parameters 
#########################################################################
go_process_of_interest <- 'pdh'
output_dir_plot <- "C:/Users/siljeew/Master_project/snRNAseq/Plots/Gene_process"

#Define outputs for t-test
output_dir_data <- "C:/Users/siljeew/Master_project/snRNAseq/Data"
output_csv_path <- "C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/"
output_xlsx_path <- "C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/"


# # Ensure output plot directory exists
# if (!dir.exists(output_dir_plot)) {
#   dir.create(output_dir_plot, recursive = TRUE)
# }
# 
# # Ensure output plot directory exists
# if (!dir.exists(output_csv_path)) {
#   dir.create(output_csv_path, recursive = TRUE)
# }
# 
# # Ensure output plot directory exists
# if (!dir.exists(output_xlsx_path)) {
#   dir.create(output_xlsx_path, recursive = TRUE)
# }
# 
# if (!dir.exists(output_dir_data)) {
#   dir.create(output_dir_data, recursive = TRUE)
# }

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




### Process seurat objects
#########################################################################
process_expression_data <- function(seurat_object, genes_of_interest) {
  data <- GetAssayData(seurat_object, slot = "data") %>% as.data.frame()
  
  # Keep only the genes from genes_of_interest
  data <- data[rownames(data) %in% toupper(genes_of_interest), ]
  rownames(data) <- toupper(rownames(data))
  
  # Convert data to long format and ensure Timepoint is included
  data_long <- data %>%
    rownames_to_column(var = "gene") %>%
    pivot_longer(cols = -gene, names_to = "cell", values_to = "expression") %>%
    mutate(
      Timepoint = seurat_object@meta.data$Timepoint[match(cell, rownames(seurat_object@meta.data))],
      Stress_Status = seurat_object@meta.data$Stress_Status[match(cell, rownames(seurat_object@meta.data))]
    )
  
  return(data_long)
}

# Process each Seurat object and return combined data
process_seurat_object <- function(file_path, genes_of_interest) {
  seurat_object <- readRDS(file_path)
  data_long <- process_expression_data(seurat_object, genes_of_interest)
  return(data_long)
}

# Process and combine data from all Seurat objects
combined_data_long <- do.call(rbind, lapply(file_paths, process_seurat_object, genes_of_interest = genes_of_interest))

head(combined_data_long)





# ### Perform anova and Tukey HSD
# ########################################################################
# Run the function
anova_tukey_results <- perform_anova_tukey(
  data = combined_data_long, 
  Stress_Status = "Stress_Status",  
  expression = "expression",        
  Timepoint = "Timepoint",          
  go_process_of_interest = go_process_of_interest,
  output_csv_path = "C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/anova_tukey_results.csv",
  output_xlsx_path = "C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/anova_tukey_results.xlsx"
)

# Inspect results
tail(anova_tukey_results, 6)






### Process seurat objects
#########################################################################
# Aggregate gene expressions by Timepoint and stress_status
go_process_aggregated <- combined_data_long %>%
  group_by(gene, Timepoint, Stress_Status) %>%
  summarize(mean_expression = mean(expression, na.rm = TRUE)) %>%
  ungroup()

# Summarize overall level by averaging gene expression levels
go_process_summary <- go_process_aggregated %>%
  group_by(Timepoint, Stress_Status) %>%
  summarize(mean_expression = mean(mean_expression, na.rm = TRUE), .groups = 'drop')

# Ensure 'stress_status' is correctly set as a factor
go_process_summary$Timepoint <- factor(go_process_summary$Timepoint, levels = unique(go_process_summary$Timepoint))
go_process_summary$Stress_Status <- factor(go_process_summary$Stress_Status, levels = c("SHAM - CM", "Not stressed CM", "Stressed CM"))

View(go_process_summary)



### Generate and save process plot
#########################################################################
# Function to generate plots with significance stars
generate_plot <- function(summary_data, go_process_of_interest) {
  
  summary_data$Timepoint <- factor(summary_data$Timepoint,
                                   levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
  
  
  go_plot <- ggplot(summary_data, aes(x = Timepoint, y = mean_expression, fill = Stress_Status)) +
    geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
    labs(x = NULL, y = paste(go_process_of_interest, "Level")) +
    theme_minimal() +
    scale_fill_manual(values = c("SHAM - CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
    scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) 
  
  return(go_plot)
}


go_plot <- generate_plot(go_process_summary, go_process_of_interest)
print(go_plot)

# Save
ggsave(file.path(output_dir_plot, paste(go_process_of_interest, "_by_stress_status.png", sep = "")), plot = go_plot, width = 8, height = 6)
ggsave(file.path(output_dir_plot, paste(go_process_of_interest, "_by_stress_status.pdf", sep = "")), plot = go_plot, width = 8, height = 6)




# ### Generate and save process plot with p-values
# #########################################################################
# go_plot <- generate_plot(go_process_summary, p_values_for_plot, go_process_of_interest)
# print(go_plot)




