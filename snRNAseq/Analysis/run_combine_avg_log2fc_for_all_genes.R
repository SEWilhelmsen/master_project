# Run combine_avg_log2fc_for_all_genes.R
# Silje Wilhelmsen

# Source necessary scripts 
source("C:/Users/siljeew/snRNAseq/Analysis/load_combine_avg_log2fc_libraries.R")
source("C:/Users/siljeew/snRNAseq/Analysis/combine_avg_log2fc_for_all_genes.R")

# Combine data into one file of average log 2 fold change values for all time points
####################################################################################

# Define the path to the data files
data_path <- "C:/Users/siljeew/snRNAseq/Data/"

# List of file names and corresponding time points
file_info <- list(
  "mouse_6h_vcm_markers_all_genes.xlsx"    = '6_hours',
  "mouse_12h_vcm_markers_all_genes.xlsx"   = '12_hours',
  "mouse_1d_vcm_markers_all_genes.xlsx"    = '1_day',
  "mouse_3d_vcm_markers_all_genes.xlsx"    = '3_days',
  "mouse_1w_vcm_markers_all_genes.xlsx"    = '1_week',
  "mouse_3w_vcm_markers_all_genes.xlsx"    = '3_weeks'
)

# Process all files and save into a list
mouse_data_list <- lapply(names(file_info), function(file) {
  process_file(file, file_info[[file]])
})

# Merge all data frames on "gene" column
merged_data <- merge_data_frames(mouse_data_list)

# # Save the merged data to a CSV file
# setwd("C:/Users/siljeew/snRNAseq/Data/")
# write.csv2(merged_data, "mouse_vcm_all_genes_avg_log2fc.csv", row.names = TRUE)
# 
# # Save the merged data to an RDS file for efficient storage and retrieval in R
# saveRDS(merged_data, "mouse_vcm_all_genes_avg_log2fc.Rds")

# Convert the "gene" column to row names
merged_data <- merged_data %>% column_to_rownames(var = "gene")

View(merged_data)

# Save the merged data to a CSV file
setwd("C:/Users/siljeew/snRNAseq/Data/")
write.csv2(merged_data, "mouse_vcm_all_genes_avg_log2fc.csv", row.names = TRUE)

# Save the merged data to an RDS file for efficient storage and retrieval in R
saveRDS(merged_data, "mouse_vcm_all_genes_avg_log2fc.Rds")
