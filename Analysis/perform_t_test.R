### Student's t-test on variable_of_interest
# Input: processed_data
# Output: .csv of t-test results. Comparing mean of two groups on each time_point
## Silje Wilhelmsen

# # Load libraries
# library(dplyr)
# library(broom)
# library(writexl)

# Function to perform t-test and save results
perform_t_test <- function(processed_data, variable_of_interest, output_csv_path, output_xlsx_path) {
  
  # Split the data by time_point
  split_data <- split(processed_data, processed_data$time_point)
  
  # Initialize a dataframe to store results
  results_df <- data.frame()
  
  # Loop through each time point and perform the t-test
  for (time_point in names(split_data)) {
    data_subset <- split_data[[time_point]]
    
    variable_of_interest_sham <- data_subset %>%
      filter(condition == 'sham') %>%
      pull(!!sym(variable_of_interest))
    
    variable_of_interest_ab <- data_subset %>%
      filter(condition == 'AB') %>%
      pull(!!sym(variable_of_interest))
    
    if (length(variable_of_interest_sham) > 0 & length(variable_of_interest_ab) > 0) {
      # Perform the t-test
      t_test_result <- t.test(variable_of_interest_sham, variable_of_interest_ab)
      
      # Extract the relevant values
      t_statistic <- t_test_result$statistic
      p_value <- t_test_result$p.value
      conf_low <- t_test_result$conf.int[1]
      conf_high <- t_test_result$conf.int[2]
      mean_group1 <- mean(variable_of_interest_sham, na.rm = TRUE)
      mean_group2 <- mean(variable_of_interest_ab, na.rm = TRUE)
      
      # Create a result data frame
      test_result <- data.frame(
        variable = variable_of_interest,
        time_point = time_point,
        t_statistic = t_statistic,
        p_value = p_value,
        conf_low = conf_low,
        conf_high = conf_high,
        mean_sham = mean_group1,
        mean_ab = mean_group2, 
        stringsAsFactors = FALSE
      )
      
      # Append the results to the results data frame
      results_df <- bind_rows(results_df, test_result)
    }
  }
  
  View(results_df)
  
  # Check if the CSV file already exists
  if (file.exists(output_csv_path)) {
    existing_results <- read.csv(output_csv_path, stringsAsFactors = FALSE)
    
    # Merge new results with existing results, updating duplicates
    merged_results <- existing_results %>%
      filter(!((time_point %in% results_df$time_point) & (variable %in% results_df$variable))) %>%
      bind_rows(results_df)
  } else {
    # If no existing file, use the new results
    merged_results <- results_df
  }
  # Write the results to a CSV file
  write.csv(merged_results, output_csv_path, row.names = FALSE)
  write_xlsx(merged_results, output_xlsx_path)
  
  view(merged_results)
  
}
