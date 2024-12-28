## Fucntions to combine all genes from .xlsx, only the column:  avg_log2FC, at all time points 
# Silje Wilhelmsen


# # Load necessary libraries
# library(readxl)
# library(dplyr)
# library(tidyr)
library(xlsx)
install.packages("xlsx")


# Function to read and process each file
process_file <- function(file, time_point) {
  data <- read_excel(file.path(data_path, file))
  data <- data %>%
    mutate(gene = toupper(gene)) %>%
    dplyr::select(gene, avg_log2FC) %>%
    rename_at(vars(avg_log2FC), ~paste0("avg_log2fc_", time_point))
  return(data)
}


# Function to merge data frames
merge_data_frames <- function(data_list) {
  merged_data <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), data_list)
  return(merged_data)
}
