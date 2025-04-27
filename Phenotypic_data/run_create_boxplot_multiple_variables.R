# Create multiple bar plots

# Boxplot
# Input: dataset of small size (e.g. n=10 per group)
# Output: 
# Silje Wilhelmsen

# Preparation
###############################################################################

# Source necessary files
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/extract_and_clean_data_from_animal_overview.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/create_boxplot_multiple_variables.R")

# Load libraries
load_libraries()

# File path to the Excel file
file_path <- "C:/Users/siljeew/Master_project/snRNAseq/Animal_overview_TP.xlsx"

# Prepare data from original excel file
processed_data <- process_excel_data(file_path, output_view = TRUE) # run the function in the source
View(processed_data) # Check the results

# Define variables of interest and output directory

variable_of_interest_list <- list("LW_(mg)", "lvw_bw", "TL_(mm)", "hw_bw")

print(colnames(processed_data))

#variable_of_interest <- "LW_(mg)"  # Change this to the desired variable
output_dir <- "C:/Users/siljeew/Master_project/snRNAseq/Data"
output_dir_plot <- "C:/Users/siljeew/Master_project/snRNAseq/Plots"

# Ensure the directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create boxplot and save
###############################################################################

# Add highlighted samples
highlighted_samples <- processed_data %>%
  filter(Used_already_for == "Multiome seq") %>%
  pull(sample_id)

print(highlighted_samples)

processed_data_highlighted <- processed_data %>%
  mutate(is_highlighted = ifelse(sample_id %in% highlighted_samples, TRUE, FALSE))

# Check the values of the samples
print(n=36, filter(processed_data_highlighted, Used_already_for == "Multiome seq"))
nrow(filter(processed_data_highlighted, Used_already_for == "Multiome seq"))

print(filter(processed_data_highlighted, is_highlighted == "TRUE"))
nrow(filter(processed_data_highlighted, is_highlighted == "TRUE"))


# Combine multiple plots 
# Prepare individual plots for each variable of interest
individual_plots <- lapply(variable_of_interest_list, function(var) {
  create_boxplot(processed_data_highlighted, "condition", var)
})

# Combine the plots into one figure
combined_plot <- grid.arrange(grobs = individual_plots, ncol = 1, nrow = 4)

# Save the combined plot
ggsave(file.path(output_dir_plot, "combined_boxplot.png"), plot = combined_plot, width = 16, height = 12)