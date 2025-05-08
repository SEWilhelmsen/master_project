# Combine average two fold change values for all timepoints into one file. 
# Silje Wilhelmsen

library(readxl)
library(dplyr)
library(purrr)


# Combine data into one file for all time points by stress
####################################################################################
# Define the path to the data files
setwd("C:/Users/siljeew/Master_project/snRNAseq/")  # Replace with your directory path
data_path <- "C:/Users/siljeew/Master_project/snRNAseq/Data/"

# List of file names and corresponding time points
file_info <- list(
  "Data/markers_6h_vcm_stress.xlsx",
  "Data/markers_12h_vcm_stress.xlsx",
  "Data/markers_1d_vcm_stress.xlsx",
  "Data/markers_3d_vcm_stress.xlsx",
  "Data/markers_1w_vcm_stress.xlsx",
  "Data/markers_3w_vcm_stress.xlsx"
)


# Combine all excel files into one
combined_data <- file_info %>%
  map_dfr(~ read_excel(.x))

combined_data <- combined_data %>%
  dplyr::rename("gene" = feature)

View(combined_data)

write_xlsx(combined_data, path = "markers_data_by_stress.xlsx")


# Combine data into one file for all time points by condition
####################################################################################
# Define the path to the data files
setwd("C:/Users/siljeew/Master_project/snRNAseq/")  # Replace with your directory path
data_path <- "C:/Users/siljeew/Master_project/snRNAseq/Data/"

# List of file names and corresponding time points
file_info <- list(
  "Data/markers_6h_vcm_group.xlsx",
  "Data/markers_12h_vcm_group.xlsx",
  "Data/markers_1d_vcm_group.xlsx",
  "Data/markers_3d_vcm_group.xlsx",
  "Data/markers_1w_vcm_group.xlsx",
  "Data/markers_3w_vcm_group.xlsx"
)


# Combine all excel files into one
combined_data <- file_info %>%
  map_dfr(~ read_excel(.x))

combined_data <- combined_data %>%
  dplyr::rename("gene" = feature)

View(combined_data)

write_xlsx(combined_data, path = "markers_data_by_group.xlsx")



# Create data frame with mean expression and percentage to SHAM
####################################################################
gene_of_interest <- "CPT2"

data_for_plot <- combined_data %>%
  filter(gene == gene_of_interest) %>%
  select(gene, group, Timepoint, avgExpr)
View(data_for_plot)

data_for_plot <- data_for_plot %>%
  group_by(Timepoint) %>%
  mutate(sham_avgExpr = avgExpr[group == "SHAM - CM"],
         percentage_to_sham = (avgExpr / sham_avgExpr) * 100) %>%
  ungroup()
View(data_for_plot)

data_for_plot <- data_for_plot %>%
  rename(Stress_Status = group)

data_for_plot <- data_for_plot %>%
  mutate(Stress_Status = str_replace(Stress_Status, "SHAM - CM", "SHAM CM"))


