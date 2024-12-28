# Boxplot
# Input: dataset of small size (e.g. n=10 per group)
# Output: 
# Silje Wilhelmsen

# Preparation
###############################################################################

# Source necessary files
source("C:/Users/siljeew/snRNAseq/Analysis/extract_and_clean_data_from_animal_overview.R")
source("C:/Users/siljeew/snRNAseq/Visualization/create_boxplot_variable.R")

# Load libraries
load_libraries()

# File path to the Excel file
file_path <- "C:/Users/siljeew/snRNAseq/Animal_overview_TP.xlsx"

# Prepare data from original excel file
processed_data <- process_excel_data(file_path, output_view = TRUE) # run the function in the source
View(processed_data) # Check the results

# Define variables of interest and output directory
variable_of_interest <- "LW_(mg)"  # Change this to the desired variable
output_dir <- "C:/Users/siljeew/snRNAseq/Data"
output_dir_plot <- "C:/Users/siljeew/snRNAseq/Plots"

# Ensure the directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create boxplot
###############################################################################

# Add highligthed samples
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

# Create and save plot
boxplot <- create_boxplot(processed_data_highlighted, "condition", variable_of_interest)
ggsave(file.path(output_dir_plot, paste(variable_of_interest, "_boxplot.png", sep = "")), plot = boxplot, width = 8, height = 6)
