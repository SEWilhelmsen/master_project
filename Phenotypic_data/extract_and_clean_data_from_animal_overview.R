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
      time_point = `Duration`,
      condition = `Op status`,
      bw = `Weight Op (g)`,
      hw = `HW (mg)`,
      hw_bw = `HW/BW`,
      rvw_bw = `RVW/BW`,
      lvw_bw = `LVW/BW`,
      lvw_hw = 'LVW/HW'
    ) %>%
    rename_with(~ gsub(" ", "_", .)) %>%
    mutate(condition = replace(condition, condition == "ORAB 0.66", "ORAB"),
           condition = replace(condition, condition == "sham", "SHAM")) %>%
    mutate(time_point = case_when(
      time_point == "1w" ~ "1 Week",
      time_point == "3w" ~ "3 Weeks",
      time_point == "3d" ~ "3 Days",
      time_point == "1d" ~ "1 Day",
      time_point == "6h" ~ "6 Hours",
      time_point == "12h" ~ "12 Hours",
      TRUE ~ time_point  # Preserve other values as they are
    )) %>%
    mutate(time_point = factor(time_point, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")))%>%
    filter(!is.na(time_point))  # Filter out rows with NA in time_point
  
  return(cleaned_data)
  
}


# Master function to process the excel data
process_excel_data <- function(file_path, sheet = 1, output_view = FALSE) {
  # Columns to extract
  columns_to_extract <- c("#Animal", "Duration", "Op status", "Weight Op (g)", "HW (mg)", "RVW (mg)", "LVW (mg)", "LW (mg)", "TL (mm)", "HW/BW", "RVW/BW", "LVW/BW", "Intention", "Used already for", "LVW/HW")
  
  # Step 1: Read the Excel file
  data <- read_excel_file(file_path, sheet)
  
  # Step 2: Extract relevant columns
  extracted_data <- extract_columns(data, columns_to_extract)
  
  # Step 3: Filter out rows with NA values in the specific column
  filtered_data <- filter_data(extracted_data, "LVW (mg)")
  
  # Step 4: Rename and clean columns
  cleaned_data <- rename_and_clean_columns(filtered_data)
  
  return(cleaned_data)
}
