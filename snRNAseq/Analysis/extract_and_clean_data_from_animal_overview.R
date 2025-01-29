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
      lvw_bw = `LVW/BW`
    ) %>%
    rename_with(~ gsub(" ", "_", .)) %>%
    mutate(condition = replace(condition, condition == "ORAB 0.66", "ORAB"),
           condition = replace(condition, condition == "sham", "SHAM")) %>%
    mutate(time_point = case_when(
      time_point == "1w" ~ "1_week",
      time_point == "3w" ~ "3_weeks",
      time_point == "3d" ~ "3_days",
      time_point == "1d" ~ "1_day",
      time_point == "6h" ~ "6_hours",
      time_point == "12h" ~ "12_hours",
      TRUE ~ time_point  # Preserve other values as they are
    )) %>%
    mutate(time_point = factor(time_point, levels = c("6_hours", "12_hours", "1_day", "3_days", "1_week", "3_weeks")))%>%
    filter(!is.na(time_point))  # Filter out rows with NA in time_point
  
  return(cleaned_data)
  
}


# Master function to process the excel data
process_excel_data <- function(file_path, sheet = 1, output_view = FALSE) {
  # Columns to extract
  columns_to_extract <- c("#Animal", "Duration", "Op status", "Weight Op (g)", "HW (mg)", "RVW (mg)", "LVW (mg)", "LW (mg)", "TL (mm)", "HW/BW", "RVW/BW", "LVW/BW", "Intention", "Used already for")
  
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
