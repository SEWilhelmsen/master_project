# Combine average two fold change values for all timepoints into one file. 
# Silje Wilhelmsen


# Source necessary scripts 
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/load_libraries_combine_avg_log2fc.R")
# source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/combine_avg_log2fc_for_all_genes.R")



# Combine data into one file of average log 2 fold change values for all time points
####################################################################################
# Define the path to the data files
data_path <- "C:/Users/siljeew/Master_project/snRNAseq/Data/"

# List of file names and corresponding time points
file_info <- list(
  "markers_6h_vcm_group.xlsx"    = '6_hours',
  "markers_12h_vcm_group.xlsx"   = '12_hours',
  "markers_1d_vcm_group.xlsx"    = '1_day',
  "markers_3d_vcm_group.xlsx"    = '3_days',
  "markers_1w_vcm_group.xlsx"    = '1_week',
  "markers_3w_vcm_group.xlsx"    = '3_weeks'
)


# Function to read and process each file
process_file <- function(file, time_point) {
  data <- read_excel(file.path(data_path, file))
  data <- data %>%
    filter(group == "ORAB") %>%
    mutate(gene = toupper(feature)) %>%
    dplyr::select(gene, logFC) %>%
    rename_at(vars(logFC), ~paste0("logFC", time_point))
  return(data)
}


# Function to merge data frames
merge_data_frames <- function(data_list) {
  merged_data <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), data_list)
  return(merged_data)
}



# Process all files and save into a list
mouse_data_list <- lapply(names(file_info), function(file) {
  process_file(file, file_info[[file]])
})

# Merge all data frames on "gene" column
merged_data <- merge_data_frames(mouse_data_list)

View(merged_data)

# Convert the "gene" column to row names
merged_data <- merged_data %>% column_to_rownames(var = "gene")

View(merged_data)

# Save the merged data to a CSV file
write.csv2(merged_data, file = file.path(data_path, "mouse_vcm_all_genes_avg_logFC.csv"), row.names = TRUE)


# Save the merged data to an RDS file for efficient storage and retrieval in R
saveRDS(merged_data, file = file.path(data_path, "mouse_vcm_all_genes_avg_logFC.Rds"))
