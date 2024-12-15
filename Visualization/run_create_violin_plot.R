# Run create violin plot to investigate variance of a gene in the different groups and at different time points
# Silje Wilhelmsen

# Preparation
##############################################################################

# Source necessary files
source("C:/Users/siljeew/snRNAseq/Visualization/load_libraries_run_create_violin_plot.R")
source("C:/Users/siljeew/snRNAseq/Analysis/process_seurat_data.R")
source("C:/Users/siljeew/snRNAseq/Visualization/create_violin_plot.R")


# Define parameters
genes_of_interest <- c('GCK')

output_dir_plot <- "C:/Users/siljeew/snRNAseq/Plots"

# Ensure output plot directory exists
if (!dir.exists(output_dir_plot)) {
  dir.create(output_dir_plot, recursive = TRUE)
}

# Define file paths for Seurat objects
file_paths <- list(
  "C:/Users/siljeew/snRNAseq/tmp/mouse_vcm_3d.Rds",
  "C:/Users/siljeew/snRNAseq/tmp/mouse_vcm_1w.Rds",
  "C:/Users/siljeew/snRNAseq/tmp/mouse_vcm_3w.Rds"
)



# Prepare data
#############################################################################

# Process and combine data from all Seurat objects
combined_data_long <- do.call(rbind, lapply(file_paths, process_seurat_object, genes_of_interest = genes_of_interest))


# Ensure 'condition' is correctly set as a factor
combined_data_summary$time_point <- factor(combined_data_summary$time_point, levels = unique(combined_data_summary$time_point))
combined_data_summary$condition <- factor(combined_data_summary$condition, levels = c("SHAM", "AB"))

# Load libraries
load_libraries()


# Create violin plot
##############################################################################
violin_plot <- create_violin_plot(combined_data_long, "condition", genes_of_interest[1])
print(violin_plot)
