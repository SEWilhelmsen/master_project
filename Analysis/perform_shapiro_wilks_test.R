### Shapiro-Wilk's test
# Input: processed_data and variable_of_interest
# Output: 
# The Shapiro-Wilk's test is used to test for normality in one sample (not comparison between several samples).  
# H0: sample distribution is normal. 
## Silje Wilhelmsen

# # Load Required Libraries
# library(dplyr)
# library(readxl)

# Function to perform Shapiro-Wilk test and save results
perform_shapiro_wilks_test <- function(processed_data, variable_of_interest, output_csv_path) {
  # Extract unique time points from the dataset
  time_points <- unique(processed_data$time_point)
  
  # Initialize a list to store the results
  results_list <- list()
  
  # Loop through each time point and perform the Shapiro-Wilk test
  for (time_point in time_points) {
    # Separate the groups by condition at the current time point
    variable_of_interest_sham <- filter(processed_data, condition == 'sham', time_point == time_point)[[variable_of_interest]]
    variable_of_interest_ab <- filter(processed_data, condition == 'AB', time_point == time_point)[[variable_of_interest]]
    
    # Ensure the extracted data are numeric and non-empty
    variable_of_interest_sham <- as.numeric(variable_of_interest_sham)
    variable_of_interest_ab <- as.numeric(variable_of_interest_ab)
    
    if (length(variable_of_interest_sham) > 0 & length(variable_of_interest_ab) > 0) {
      # Perform the Shapiro-Wilk test for sham group
      sham_shapiro_wilks_result <- shapiro.test(variable_of_interest_sham)
      # Record results for sham group
      sham_result <- data.frame(
        time_point = time_point,
        group = 'sham',
        w_statistic = sham_shapiro_wilks_result$statistic,
        p_value = sham_shapiro_wilks_result$p.value,
        mean = mean(variable_of_interest_sham),
        variable_of_interest = variable_of_interest
      )
      
      # Perform the Shapiro-Wilk test for ab group
      ab_shapiro_wilks_result <- shapiro.test(variable_of_interest_ab)
      # Record results for ab group
      ab_result <- data.frame(
        time_point = time_point,
        group = 'AB',
        w_statistic = ab_shapiro_wilks_result$statistic,
        p_value = ab_shapiro_wilks_result$p.value,
        mean = mean(variable_of_interest_ab),
        variable_of_interest = variable_of_interest
      )
      
      # Append results to the list
      results_list <- append(results_list, list(sham_result, ab_result))
    }
  }   
  
  # Convert results list to data frame
  shapiro_wilks_results_df <- do.call(rbind, results_list)
  
  # Print the results data frame
  print(shapiro_wilks_results_df)
  
  # Check the results
  View(shapiro_wilks_results_df)
  
  # Load existing results if the file exists
  if (file.exists(output_csv_path)) {
    existing_results <- read.csv(output_csv_path)
    
    # Merge new results with existing results, updating duplicates
    merged_results <- existing_results %>%
      filter(!(time_point %in% shapiro_wilks_results_df$time_point & 
                 group %in% shapiro_wilks_results_df$group & 
                 variable_of_interest %in% shapiro_wilks_results_df$variable_of_interest)) %>%
      bind_rows(shapiro_wilks_results_df)
  } else {
    # If no existing file, use the new results
    merged_results <- shapiro_wilks_results_df
  }
  
  # Write the merged results to CSV file
  write.csv(merged_results, output_csv_path, row.names = FALSE)
}

