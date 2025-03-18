# Perform ANOVA and Tukey HSD on Stress_status groups 
# Comparison of RNA seq for gene ontology process
# Silje Wilhelmsen

# Function to perform anova and Tukey
perform_anova_tukey <- function(data, Stress_Status, expression, Timepoint, go_process_of_interest, output_csv_path, output_xlsx_path) {
  # Ensure column names exist and are strings
  stopifnot(is.character(Stress_Status), is.character(expression), is.character(Timepoint))
  stopifnot(all(c(Stress_Status, expression, Timepoint) %in% names(data)))
  
  # Split data by the given timepoint
  split_data <- split(data, data[[Timepoint]])
  results <- list()
  
  # Conduct ANOVA and Tukey's test
  for (time_point in names(split_data)) {
    data_subset <- split_data[[time_point]]
    
    if (nrow(data_subset) >= 3) {  # Ensure enough data for ANOVA
      # Perform ANOVA
      anova_result <- aov(expression ~ Stress_Status, data = data_subset)
      
      # Tukey HSD post-hoc testing
      tukey_result <- TukeyHSD(anova_result)
      
      # Access the correct matrix from Tukey result
      tukey_df <- as.data.frame(tukey_result[["Stress_Status"]])
      
      # Add rownames as a column for group comparisons
      tukey_df$comparison <- rownames(tukey_df)
      
      # Add process and timepoint to results
      tukey_df$process <- go_process_of_interest
      tukey_df$Timepoint <- time_point
      
      # Collect results
      results[[time_point]] <- tukey_df
    }
  }
  
  results_df <- bind_rows(results)
  
  # Update existing results
  if (file.exists(output_csv_path)) {
    existing_results <- read_csv(output_csv_path, show_col_types = FALSE)
    
    # Remove old entries for the current process and time points, then merge
    updated_results <- existing_results %>%
      filter(!(process == go_process_of_interest & Timepoint %in% results_df$Timepoint))
    
    results_df <- bind_rows(updated_results, results_df)
  }
  
  # Write updated results
  write_csv(results_df, output_csv_path)
  write.xlsx(results_df, output_xlsx_path)
  
  return(results_df)
}