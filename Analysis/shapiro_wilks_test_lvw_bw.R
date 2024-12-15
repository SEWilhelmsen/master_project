### Shapiro-Wilk's test 
# Input: animal_data_rnaseq_renamed, variable: LVW/BW
# Output: 
# The Shapiro-Wilk's test is used to test for normality in one sample (not comparison between several samples).  
# H0: sample distribution is normal. 
## Silje Wilhelmsen

# Load Required Libraries
load_libraries <- function() {
  library(readxl)
  library(dplyr)
}

# Function to perform Shapiro-Wilk test and save results
perform_shapiro_wilks_test <- function(data, variable, output_csv_path) {
  # Extract unique time points from the dataset
  time_points <- unique(data$time_point)
  
  # Initialize a list to store the results
  results_list <- list()
  
  # Loop through each time point and perform the Shapiro-Wilk test
  for (time_point in time_points) {
    # Separate the groups by condition at the current time point
    variable_sham <- filter(data, condition == 'sham', time_point == time_point)[[variable]]
    variable_ab <- filter(data, condition == 'AB', time_point == time_point)[[variable]]
    
    # Perform the Shapiro-Wilk test for sham group
    sham_shapiro_wilks_result <- shapiro.test(variable_sham)
    
    # Record results for sham group
    sham_result <- data.frame(
      time_point = time_point,
      group = 'sham',
      w_statistic = sham_shapiro_wilks_result$statistic,
      p_value = sham_shapiro_wilks_result$p.value,
      mean = mean(variable_sham)
    )
    
    # Perform the Shapiro-Wilk test for ab group
    ab_shapiro_wilks_result <- shapiro.test(variable_ab)
    
    # Record results for ab group
    ab_result <- data.frame(
      time_point = time_point,
      group = 'ab',
      w_statistic = ab_shapiro_wilks_result$statistic,
      p_value = ab_shapiro_wilks_result$p.value,
      mean = mean(variable_ab)
    )
    
    # Append results to the list
    results_list <- append(results_list, list(sham_result, ab_result))
  }
  
  # Convert results list to data frame
  shapiro_wilks_results_df <- do.call(rbind, results_list)
  
  # Print the results data frame
  print(shapiro_wilks_results_df)
  
  # # Write the results to a CSV file
  # write.csv(shapiro_wilks_results_df, output_csv_path, row.names = FALSE)
}

View(processed_data)


filtered_data <- processed_data %>%
  filter(condition == "sham" & time_point == "3_weeks") %>%
  select(sample_id, condition, time_point, lvw_bw)
View(filtered_data)

shapiro.test(filtered_data$lvw_bw)
