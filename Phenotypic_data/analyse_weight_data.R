# Create barplot of weight data

# Load libraries 
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyverse)
library(viridis)
library(hrbrthemes)

# Import data
Animal_overview_TP_copy <- read_excel("C:/Users/siljeew/Master_project/Phenotypic_data/Animal_overview_TP.xlsx")


# Extract relevant data 
animal_data <- Animal_overview_TP_copy[, c("#Animal", "Duration", "Op status", "HW (mg)", "RVW (mg)", "LVW (mg)", "LW (mg)", "TL (mm)", "HW/BW", "RVW/BW","LVW/BW", "Intention", "Used already for", "LVW/HW")]
#View(animal_data)

## Filter for value in a column
#filtered_df <- subset(df, Used.for == "Multiome seq")
# 
# # Remove samples with missing data (NA) in the LVW (mg) column
# animal_data_rnaseq_filtered <- na.omit(animal_data_rnaseq)

# Remove rows with NA values in specific columns
animal_data_filtered <- animal_data %>% filter(!is.na(`LVW (mg)`))

# Change naming to match other scripts
animal_data_renamed <- animal_data_filtered %>%
  # Rename specific columns
  rename(
    time_point = `Duration`,
    condition = `Op status`,
    sample_id = `#Animal`
  ) %>%
  # Remove spaces from all column names
  rename_with(~ gsub(" ", "_", .)) %>%
  mutate(condition = replace(condition, condition == "ORAB 0.66", "ORAB"), 
         condition = replace(condition, condition == "sham", "SHAM")) %>%
  
  # Replace specific values in the time_point column using case_when
  mutate(time_point = case_when(
    time_point == "1w" ~ "1 Week",
    time_point == "3w" ~ "3 Weeks",
    time_point == "3d" ~ "3 Days",
    time_point == "1d" ~ "1 Day",
    time_point == "6h" ~ "6 Hours",
    time_point == "12h" ~ "12 Hours",
    TRUE ~ time_point  # Preserve other values as they are
  )) %>%
  # Convert time_point to a factor and set the desired order of its levels
  mutate(time_point = factor(time_point, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")))

# View the updated data frame
View(animal_data_renamed)

# Define samples for highlighting
highlighted_samples <- animal_data_renamed %>%
  filter(Used_already_for == "Multiome seq") %>%
  pull(sample_id)


