# Kruskal wallis test 
# Non parametric test of multiple groups. 
# H0: the groups are equal 

# library(dplyr)

# Function to perform Kruskal-Wallis test and save the results
perform_kruskal_wallis_test <- function(processed_data, variable, output_path) {
  
  # Create a new column 'group' by combining 'condition' and 'time_point'
  processed_data <- processed_data %>%
    mutate(group = paste(condition, time_point, sep = "_"))
  
  result <- kruskal.test(as.formula(paste(variable, "~ group")), data = processed_data)
  
  # Convert the htest object to a data frame
  result_df <- data.frame(
    variable_of_interest = variable,
    statistic = result$statistic,
    parameter = result$parameter,
    p_value = result$p.value,
    method = result$method,
    data_name = result$data.name,
    stringsAsFactors = FALSE
  )
  
  View(result_df)
  
  # Check if the CSV file already exists
  if (file.exists(output_path)) {
    existing_results <- read.csv(output_path)
    
    # Merge new results with existing results, updating duplicates
    merged_results <- existing_results %>%
      filter(!(variable_of_interest %in% result_df$variable_of_interest)) %>%
      bind_rows(result_df)
  } else {
    # If no existing file, use the new results
    merged_results <- result_df
  }
  
  # Save the results to a CSV file
  write.csv(merged_results, output_path, row.names = FALSE)
  
  View(merged_results)
  
  # Optional: View the results
  return(merged_results)
}
