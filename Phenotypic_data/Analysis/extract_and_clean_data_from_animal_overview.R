# Extract relevant data and prepare for hypothesis testing
# Silje Wilhelmsen

# Load libraries
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyverse)
library(viridis)
library(hrbrthemes)

# Function to read an Excel file. Remember to define file_path in your running script
read_excel_file <- function(file_path, sheet = 1) {
  data <- read_excel(file_path, sheet = sheet)
  return(data)
}

# Function to extract specific columns
extract_columns <- function(data, columns) {
  extracted_data <- data %>% select(all_of(columns))
  return(extracted_data)
}

# Function to filter data based on NA values in a column
filter_data <- function(data, column) {
  filtered_data <- data %>% filter(!is.na(.data[[column]]))
  return(filtered_data)
}

# Function to rename and clean columns
rename_and_clean_columns <- function(data) {
  cleaned_data <- data %>%
    rename(
      sample_id = `#Animal`,
      Timepoint = `Duration`,
      Condition = `Op status`,
      bw = `Weight Op (g)`,
      hw = `HW (mg)`,
      hw_bw = `HW/BW`,
      rvw_bw = `RVW/BW`,
      lvw_bw = `LVW/BW`,
      lvw_hw = 'LVW/HW'
    ) %>%
    rename_with(~ gsub(" ", "_", .)) %>%
    mutate(Condition = replace(Condition, Condition == "ORAB 0.66", "ORAB"),
           Condition = replace(Condition, Condition == "sham", "SHAM")) %>%
    mutate(Timepoint = case_when(
      Timepoint == "1w" ~ "1 Week",
      Timepoint == "3w" ~ "3 Weeks",
      Timepoint == "3d" ~ "3 Days",
      Timepoint == "1d" ~ "1 Day",
      Timepoint == "6h" ~ "6 Hours",
      Timepoint == "12h" ~ "12 Hours",
      TRUE ~ Timepoint  # Preserve other values as they are
    )) %>%
    mutate(Timepoint = factor(Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")))%>%
    filter(!is.na(Timepoint))  # Filter out rows with NA in Timepoint
  
  return(cleaned_data)
  
}


# Master function to process the excel data
process_excel_data <- function(file_path, sheet = 1, output_view = FALSE) {
  columns_to_extract <- c("#Animal", "Duration", "Op status", "Weight Op (g)", "HW (mg)", "RVW (mg)", "LVW (mg)", "LW (mg)", "TL (mm)", "HW/BW", "RVW/BW", "LVW/BW", "Intention", "Used already for", "LVW/HW")
  
  # Read the Excel file
  data <- read_excel_file(file_path, sheet)
  
  # Extract relevant columns
  extracted_data <- extract_columns(data, columns_to_extract)
  
  # Filter out rows with NA values in the specific column
  filtered_data <- filter_data(extracted_data, "LVW (mg)")
  
  # Rename and clean columns
  cleaned_data <- rename_and_clean_columns(filtered_data)
  
  return(cleaned_data)
}
