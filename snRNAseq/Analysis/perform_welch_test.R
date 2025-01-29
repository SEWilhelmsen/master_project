# Welch test for unequal variance between the groups 
# Input: processed_data
# Output: .csv of t-test results. Comparing ...
## Silje Wilhelmsen

# # Load libraries
# library(dplyr)
# library(broom)

# Function to filter data and perform t-test
# Function to filter data and perform Welch's t-test for a given time point
# perform_welch_test_for_time_point <- function(data_subset, variable_of_interest) {
#   variable_of_interest_sham <- data_subset %>%
#     filter(condition == 'sham') %>%
#     pull(!!sym(variable_of_interest))
#   
#   variable_of_interest_ab <- data_subset %>%
#     filter(condition == 'AB') %>%
#     pull(!!sym(variable_of_interest))
#   
#   if (length(variable_of_interest_sham) > 0 & length(variable_of_interest_ab) > 0) {
#     # Perform the Welch's t-test
#     welch_test_result <- t.test(variable_of_interest_sham, variable_of_interest_ab, var.equal = FALSE)
#     
#     # Extract the relevant values into a data frame
#     return(data.frame(
#       variable = variable_of_interest,
#       time_point = unique(data_subset$time_point),
#       t_statistic = welch_test_result$statistic,
#       p_value = welch_test_result$p.value,
#       conf_low = welch_test_result$conf.int[1],
#       conf_high = welch_test_result$conf.int[2],
#       mean_sham = mean(variable_of_interest_sham, na.rm = TRUE),
#       mean_ab = mean(variable_of_interest_ab, na.rm = TRUE),
#       stringsAsFactors = FALSE
#     ))
#   } else {
#     return(NULL)
#   }
# }
# 
# # Main function to perform Welch's t-test and save results
# perform_welch_test <- function(processed_data, variable_of_interest, output_csv_path, output_xlsx_path) {
#   # Create a new column 'group' by combining 'condition' and 'time_point'
#   processed_data <- processed_data %>%
#     mutate(group = paste(condition, time_point, sep = "_"))
#   
#   # Perform the Welch's t-test for each time point
#   results_df <- split(processed_data, processed_data$time_point) %>%
#     lapply(perform_welch_test_for_time_point, variable_of_interest = variable_of_interest) %>%
#     bind_rows()
#   
#   if (nrow(results_df) > 0) {
#     # View the newly created results
#     View(results_df)
#     
#     # Check if the CSV file already exists
#     if (file.exists(output_csv_path)) {
#       existing_results <- tryCatch({
#         read.csv(output_csv_path, stringsAsFactors = FALSE)
#       }, error = function(e) {
#         # Handle errors during reading
#         message("Error reading the existing results file: ", e)
#         data.frame()
#       })
#       
#       # Merge new results with existing results, updating duplicates
#       merged_results <- existing_results %>%
#         filter(!(time_point %in% results_df$time_point & variable %in% results_df$variable)) %>%
#         bind_rows(results_df)
#     } else {
#       # If no existing file, use the new results
#       merged_results <- results_df
#     }
#     
#     # Write the merged results to CSV and Excel files
#     write.csv(merged_results, output_csv_path, row.names = FALSE)
#     write_xlsx(merged_results, output_xlsx_path)
#     
#     # Optional: View merged results at the end
#     View(merged_results)
#     
#     # Return the merged results
#     return(merged_results)
#   } else {
#     message("No valid Welch test results were generated.")
#     return(NULL)
#   }
# }




# Function to filter data and perform Welch's t-test for a given time point
perform_welch_test_for_time_point <- function(data_subset, variable_of_interest) {
  variable_of_interest_sham <- data_subset %>%
    filter(condition == 'sham') %>%
    pull(!!sym(variable_of_interest))
  
  variable_of_interest_ab <- data_subset %>%
    filter(condition == 'AB') %>%
    pull(!!sym(variable_of_interest))
  
  if (length(variable_of_interest_sham) > 0 & length(variable_of_interest_ab) > 0) {
    # Perform the Welch's t-test
    welch_test_result <- t.test(variable_of_interest_sham, variable_of_interest_ab, var.equal = FALSE)
    
    # Extract the relevant values into a data frame
    return(data.frame(
      variable = variable_of_interest,
      time_point = unique(data_subset$time_point),
      t_statistic = welch_test_result$statistic,
      p_value = welch_test_result$p.value,
      conf_low = welch_test_result$conf.int[1],
      conf_high = welch_test_result$conf.int[2],
      mean_sham = mean(variable_of_interest_sham, na.rm = TRUE),
      mean_ab = mean(variable_of_interest_ab, na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
  } else {
    return(NULL)
  }
}


# Main function to perform Welch's t-test and save results
perform_welch_test <- function(processed_data, variable_of_interest, output_csv_path, output_xlsx_path) {
  # Create a new column 'group' by combining 'condition' and 'time_point'
  processed_data <- processed_data %>%
    mutate(group = paste(condition, time_point, sep = "_"))
  
  # Perform the Welch's t-test for each time point
  results_df <- split(processed_data, processed_data$time_point) %>%
    lapply(perform_welch_test_for_time_point, variable_of_interest = variable_of_interest) %>%
    bind_rows()
  
  if (nrow(results_df) > 0) {
    # View the newly created results
    View(results_df)
    
    # Check if the CSV file already exists
    if (file.exists(output_csv_path)) {
      existing_results <- tryCatch({
        read.csv(output_csv_path, stringsAsFactors = FALSE)
      }, error = function(e) {
        # Handle errors during reading
        message("Error reading the existing results file: ", e)
        data.frame()
      })
      
      # Merge new results with existing results, updating duplicates
      merged_results <- existing_results %>%
        filter(!(time_point %in% results_df$time_point & variable %in% results_df$variable)) %>%
        bind_rows(results_df)
    } else {
      # If no existing file, use the new results
      merged_results <- results_df
    }
    
    # Write the merged results to CSV and Excel files
    write.csv(merged_results, output_csv_path, row.names = FALSE)
    write_xlsx(merged_results, output_xlsx_path)
    
    # Optional: View merged results at the end
    View(merged_results)
    
    # Return the merged results
    return(merged_results)
  } else {
    message("No valid Welch test results were generated.")
    return(NULL)
  }
}




