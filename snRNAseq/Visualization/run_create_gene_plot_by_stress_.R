# Create gene plot by stress status
# Define gene of interest, read data, perform ANOVA+TukeyHSD, save results and prepare data
# for temporal line plot.  
# Silje Wilhelmsen


source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/load_libraries_run_create_barplot_go_process.R")
# source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/define_go_terms_and_find_genes.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/perform_anova_tukey_one_gene_by_stress_status.R")

# file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/perform_anova_tukey_one_gene_by_stress_status.R")



# Specify the gene you’re interested in 
gene_of_interest <- "CPT1A"  # Change to gene of interest

# Specify paths for ANOVA and Tukey HSD results
output_csv_path <- "C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/anova_tukey_gene_results_new.csv"
output_xlsx_path <- "C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/anova_tukey_gene_results_new.xlsx"


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



### Process seurat objects
#########################################################################
process_expression_data <- function(seurat_object, gene_of_interest) {
  data <- GetAssayData(seurat_object, layer = "data") %>% as.data.frame()
  
  # Keep only the gene of interest
  data <- data[rownames(data) == toupper(gene_of_interest), , drop = FALSE]
  rownames(data) <- toupper(rownames(data)) # Make sure gene name consistency
  
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
process_seurat_object <- function(file_path, gene_of_interest) {
  seurat_object <- readRDS(file_path)
  data_long <- process_expression_data(seurat_object, gene_of_interest)
  return(data_long)
}

# Process and combine data from all Seurat objects
combined_data_long <- do.call(rbind, lapply(file_paths, process_seurat_object, gene_of_interest = gene_of_interest))
combined_data_long <- combined_data_long %>%
  mutate(Stress_Status = str_replace(Stress_Status, "SHAM - CM", "SHAM CM"))


write.csv(combined_data_long, "C:/Users/siljeew/Master_project/snRNAseq/tmp/combined_data_long_cpt1a.csv") # Change gene



# Perform anova and TukeyHSD
#########################################################################
# Run the function
anova_tukey_gene_results <- perform_anova_tukey_for_gene(
  data = combined_data_long,
  Stress_Status = "Stress_Status",
  expression = "expression",
  Timepoint = "Timepoint",
  gene_of_interest = gene_of_interest,
  output_csv_path = output_csv_path,
  output_xlsx_path = output_xlsx_path
)



# Prepare data for plot
########################################################################
# Aggregate gene expressions by timepoint and Stress status
go_process_aggregated <- combined_data_long %>%
  group_by(gene, Timepoint, Stress_Status) %>%
  summarize(mean_expression = mean(expression, na.rm = TRUE)) %>%
  ungroup()

# Summarize overall level by averaging gene expression levels
go_process_summary <- go_process_aggregated %>%
  group_by(Timepoint, gene, Stress_Status) %>%
  summarize(mean_expression = mean(mean_expression, na.rm = TRUE), .groups = 'drop')

# Ensure 'condition' is correctly set as a factor
go_process_summary$Timepoint <- factor(go_process_summary$Timepoint, levels = unique(go_process_summary$Timepoint))
go_process_summary$Stress_Status <- factor(go_process_summary$Stress_Status, levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))


# Calculate percentages relative to SHAM for each time point
go_process_summary <- go_process_summary %>%
  group_by(Timepoint) %>%
  mutate(sham_mean_expression = mean_expression[Stress_Status == "SHAM CM"],
         percentage_to_sham = mean_expression / sham_mean_expression * 100) %>%
  ungroup()
head(go_process_summary)



# ANOVA and Benjamini-Hochberg correction
# ###############################################################################
# anova_bh_gene_results <- perform_anova_bh_for_gene(
#   data = combined_data_long, 
#   Stress_Status = "Stress_Status",  
#   expression = "expression",        
#   Timepoint = "Timepoint",          
#   gene_of_interest = gene_of_interest,
#   output_csv_path = "C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/anova_bh_gene_results.csv",
#   output_xlsx_path = "C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/anova_bh_gene_results.xlsx"
# )
# 
# print(anova_bh_gene_results$pairwise_comparison)



data_for_t_test <- combined_data_long %>%
  mutate(Group_Timepoint = paste(Stress_Status, Timepoint, sep = "_"))

pairwise_results <- pairwise.t.test(data_for_t_test$expression, data_for_t_test$Group_Timepoint, p.adjust.method = "none", pool.sd = FALSE, alternative = c("two.sided"))
print(pairwise_results)

pairwise_results_bh <- pairwise.t.test(data_for_t_test$expression, data_for_t_test$Group_Timepoint, p.adjust.method = "BH", pool.sd = FALSE, alternative = c("two.sided"))
print(pairwise_results_bh)

# Manually copy the p-values into "C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/t_test_gene_by_stress.xlsx"
library(readxl)
t_test_gene_by_stress <- read_excel("Data/Stress_Status/t_test_gene_by_stress.xlsx")
# View(t_test_gene_by_stress)

# Use this for labeling
t_test_gene_by_stress <- t_test_gene_by_stress %>%
  mutate(significance_label = ifelse(group1 != "SHAM CM" & group2 != "SHAM CM" & p_adj < 0.05, "\u25B2",
                                     ifelse(p_adj < 0.001, "***",
                                            ifelse(p_adj < 0.01, "**",
                                                   ifelse(p_adj < 0.05, "*", "ns")))))

# go_process_summary_test <- go_process_summary %>%
#   left_join(
#     t_test_gene_by_stress %>% dplyr::select(gene, Timepoint, Group2, p_adj), 
#     by = c("gene" = "gene", "Timepoint" = "Timepoint", "Stress_Status" = "Group2"))
# 
# View(go_process_summary_test)



# Combine ANOVA and Tukey results with mean expression
###############################################################################
# Convert factors to characters 
anova_tukey_gene_results <- anova_tukey_gene_results %>%
  mutate(Timepoint = as.character(Timepoint),
         group1 = as.character(group1),
         group2 = as.character(group2))

go_process_summary <- go_process_summary %>%
  mutate(Timepoint = as.character(Timepoint),
         Stress_Status = as.character(Stress_Status))

# view(anova_tukey_gene_results)

source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/insert_mean_anova_tukey_result_by_stress.R")


#Add a new column for sham_mean_expression
anova_tukey_mean <- anova_tukey_mean %>%
  group_by(Timepoint, gene) %>%
  mutate(sham_mean_expression = group1_mean_expression[group1 == "SHAM CM"])


# Add a new column for significance label
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(significance_label = ifelse(group1 != "SHAM CM" & group2 != "SHAM CM" & `p adj` < 0.05, "\u25B2",
                                     ifelse(`p adj` < 0.001, "***",
                                            ifelse(`p adj` < 0.01, "**",
                                                   ifelse(`p adj` < 0.05, "*", "ns")))))

# Print to verify the labeling
#View(anova_tukey_mean)


source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/insert_significance_label_to_go_process_summary_by_stress.R")

write_xlsx(go_process_summary, path = "C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/go_process_summary_cpt1a.xlsx") # Change gene
# Remember to verify that the stars in the plot corresponds to the t_test_gene_by_stress p_adj 


#############################################################################################################
write_xlsx(anova_tukey_mean, path = "C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/anova_tukey_mean.xlsx") 
write.csv(anova_tukey_mean, file = "C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/anova_tukey_mean.csv")
#############################################################################################################


# To create a line plot, go here:
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/create_line_plot_one_gene_with_significance_by_stress.R")

