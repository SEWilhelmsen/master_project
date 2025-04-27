# Mann-Whitney U test for non-parametric test of difference between two independent groups
# Input:
# Output:
# Silje Wilhelmsen


# Mann-Whitney U test for significance between two independent groups,
library(dplyr)
library(readr)
library(writexl) 



# perform_mann_whitney_test <- function(processed_data, variable_of_interest, output_csv_path, output_xlsx_path) {
#   # Perform the Mann-Whitney U test for each time point
#   results_df <- split(processed_data, processed_data$Timepoint) %>%
#     lapply(perform_mann_whitney_test_for_timepoint, variable_of_interest = variable_of_interest) %>%
#     bind_rows()
#   
#   if (nrow(results_df) > 0) {
#     View(results_df)
#     
#     # Write the results to CSV and Excel files, assuming no previous results
#     write.csv(results_df, output_csv_path, row.names = FALSE)
#     write_xlsx(results_df, output_xlsx_path)
#     
#     View(results_df)
#     return(results_df)
#   } else {
#     message("No valid Mann-Whitney test results were generated.")
#     return(NULL)
#   }
# }
# 
# # Define output paths
# output_csv_path <- file.path(output_dir, "mann_whitney_results.csv")
# output_xlsx_path <- file.path(output_dir, "mann_whitney_results.xlsx")
# 
# # Perform tests on the 'bw' variable
# mann_whitney_results <- perform_mann_whitney_test(processed_data, "bw", output_csv_path, output_xlsx_path)


perform_mann_whitney_test_for_timepoint <- function(timepoint_data, variable_of_interest) {
  # Subset the data by Condition
  variable_of_interest_sham <- timepoint_data[[variable_of_interest]][timepoint_data$Condition == "SHAM"]
  variable_of_interest_orab <- timepoint_data[[variable_of_interest]][timepoint_data$Condition == "ORAB"]

  # Perform the Mann-Whitney U test
  test_result <- wilcox.test(
    as.formula(paste(variable_of_interest, "~ Condition")),
    data = timepoint_data,
    exact = FALSE  # Use exact = FALSE to handle ties gracefully
  )

  # Extract the relevant values into a data frame
  return(data.frame(
    variable = variable_of_interest,
    Timepoint = unique(timepoint_data$Timepoint),
    statistic = test_result$statistic,
    p_value = test_result$p.value,
    mean_sham = mean(variable_of_interest_sham, na.rm = TRUE),
    mean_orab = mean(variable_of_interest_orab, na.rm = TRUE)
  ))
}


perform_mann_whitney_test <- function(processed_data, variable_of_interest, output_csv_path, output_xlsx_path) {
  # Perform the Mann-Whitney U test for each time point
  results_df <- split(processed_data, processed_data$Timepoint) %>%
    lapply(perform_mann_whitney_test_for_timepoint, variable_of_interest = variable_of_interest) %>%
    bind_rows()
  
  if (nrow(results_df) > 0) {
    View(results_df)
    
    # Check if the CSV file already exists
    if (file.exists(output_csv_path)) {
      existing_results <- tryCatch({
        read.csv(output_csv_path, stringsAsFactors = FALSE)
      }, error = function(e) {
        message("Error reading the existing results file: ", e)
        data.frame()
      })
      
      # Ensure column compatibility
      if ("Timepoint" %in% names(existing_results) && "variable" %in% names(existing_results)) {
        merged_results <- existing_results %>%
          bind_rows(results_df)
      } else {
        message("Columns 'Timepoint' and 'variable' are not found in existing results, falling back to new results.")
        merged_results <- results_df
      }
    } else {
      merged_results <- results_df
    }
    
    # Write the merged results to CSV and Excel files
    write.csv(merged_results, output_csv_path, row.names = FALSE)
    write_xlsx(merged_results, output_xlsx_path)
    
    View(merged_results)
    return(merged_results)
  } else {
    message("No valid Mann-Whitney test results were generated.")
    return(NULL)
  }
}



