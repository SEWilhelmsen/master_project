# Create genes plot by stress status
# Silje Wilhelmsen

library(ggbreak)
library(stringr)


source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/load_libraries_run_create_barplot_go_process.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/define_go_terms_and_find_genes.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/perform_anova_tukey_one_gene_by_stress_status.R")

# glycolysis = toupper(c("HK1", "HK2", "GCK", "PFK1", "PFKM", "PK")),
# pdh = toupper(c("PDH", "PDHA1", "PDHA2", "DLAT", "DLD", "PDK1", "PDK2",  "PDHX")),
# glucose_transmembrane_transport = toupper(c("SLC2A4", "SLC27A1", "SLC27A4")),
# fatty_acid_transport = toupper(c("FABP3", "FABPH", "CD36", "FATP4", "CACP", "MCAT")),
# fao = toupper(c("CPT1B", "CPT1A", "ACADVL", "CPT2")),
# hydroxybutyrate_dehydrogenase_activity = toupper(c("BDH1", "BDH2")),
# ketone_catabolism = toupper(c("ACAT1", "OXCT1", "OXCT2", "OXCT2A", "OXCT2B")),

# Specify the gene you’re interested in
gene_of_interest <- "BDH1" 


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
  data <- GetAssayData(seurat_object, slot = "data") %>% as.data.frame()
  
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


head(combined_data_long)
print(gene_of_interest)
print(unique(combined_data_long$gene))

combined_data_long <- combined_data_long %>%
  mutate(Stress_Status = str_replace(Stress_Status, "SHAM - CM", "SHAM CM"))



# ### Perform anova and Tukey HSD
# ########################################################################
# Run the function
anova_tukey_gene_results <- perform_anova_tukey_for_gene(
  data = combined_data_long, 
  Stress_Status = "Stress_Status",  
  expression = "expression",        
  Timepoint = "Timepoint",          
  gene_of_interest = gene_of_interest,
  output_csv_path = "C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/anova_tukey_gene_results.csv",
  output_xlsx_path = "C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/anova_tukey_gene_results.xlsx"
)

# Inspect results
tail(anova_tukey_gene_results, 6)



head(go_process_summary)
tail(anova_tukey_gene_results)

# Prepare data for plot
########################################################################
# Aggregate gene expressions by timepoint and Stress status
go_process_aggregated <- combined_data_long %>%
  group_by(gene, Timepoint, Stress_Status) %>%
  summarize(mean_expression = mean(expression, na.rm = TRUE)) %>%
  ungroup()

# Summarize overall level by averaging gene expression levels
go_process_summary <- go_process_aggregated %>%
  group_by(Timepoint, Stress_Status) %>%
  summarize(mean_expression = mean(mean_expression, na.rm = TRUE), .groups = 'drop')

# Ensure 'condition' is correctly set as a factor
go_process_summary$Timepoint <- factor(go_process_summary$Timepoint, levels = unique(go_process_summary$Timepoint))
go_process_summary$Stress_Status <- factor(go_process_summary$Stress_Status, levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))

View(go_process_summary)


# Calculate percentages relative to SHAM for each time point
go_process_summary <- go_process_summary %>%
  group_by(Timepoint) %>%
  mutate(sham_mean_expression = mean_expression[Stress_Status == "SHAM CM"],
         percentage_to_sham = mean_expression / sham_mean_expression * 100) %>%
  ungroup()

# View the modified data frame
head(go_process_summary)





# Barplot of one gene 
###########################################
single_gene_plot <- ggplot(go_process_summary, aes(x = Timepoint, y = percentage_to_sham, fill = Stress_Status)) +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge()) +
  labs(x = NULL, y = "Expression compared to SHAM CM (%)", title = paste("Expression of", gene_of_interest)) +
  theme_minimal() +
  scale_fill_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_break(c(10, 60)) +
  scale_y_continuous(breaks = seq(0, 260, by = 20)) + 
  theme(axis.text.x = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(size = 10, colour = "black"))

single_gene_plot 


ggsave(file.path(output_dir_plot, paste(gene_of_interest, "_percentage_by_stress_status.png", sep = "")), plot = single_gene_plot, width = 8, height = 6)
ggsave(file.path(output_dir_plot, paste(gene_of_interest, "_percentage_by_stress_status.pdf", sep = "")), plot = single_gene_plot, width = 8, height = 6)


