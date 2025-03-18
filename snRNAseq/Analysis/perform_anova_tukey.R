# ANOVA and Tukey HSD test
# Input:
# Output: 
# Silje Wilhelmsen


library(openxlsx)
library(dplyr)

# Function to perform ANOVA and Tukey's HSD test for a given variable and save separate results
perform_anova_tukey_test <- function(processed_data, variable_of_interest, anova_csv_path, anova_xlsx_path, tukey_csv_path, tukey_xlsx_path) {
  
  # Update this formula to include the interaction term
  formula <- as.formula(paste(variable_of_interest, "~ Condition + Timepoint + Condition:Timepoint"))
  
  # Perform the ANOVA
  anova_result <- aov(formula, data = processed_data)
  anova_summary <- summary(anova_result)
  
  # Perform Tukey's HSD post-hoc test
  tukey_result <- TukeyHSD(anova_result, "Condition:Timepoint")
  
  # Extract ANOVA summary into a data frame
  anova_df <- data.frame(
    variable = variable_of_interest,
    term = rownames(anova_summary[[1]]),
    df = anova_summary[[1]]$Df,
    sum_sq = anova_summary[[1]]$"Sum Sq",
    mean_sq = anova_summary[[1]]$"Mean Sq",
    F_value = anova_summary[[1]]$"F value",
    p_value = anova_summary[[1]]$"Pr(>F)",
    stringsAsFactors = FALSE
  )
  
  # Extract Tukey's HSD results into a data frame
  tukey_df <- as.data.frame(tukey_result$`Condition:Timepoint`, stringsAsFactors = FALSE)
  tukey_df$comparison <- rownames(tukey_df)
  rownames(tukey_df) <- NULL
  tukey_df <- tukey_df %>%
    rename(
      diff = "diff",
      lwr = "lwr",
      upr = "upr",
      adjusted_p_value = "p adj"
    ) %>%
    mutate(variable = variable_of_interest)
  
  # Reorder columns so that 'comparison' appears first
  tukey_df <- tukey_df %>%
    select(comparison, everything())
  
  # Save the ANOVA results
  if (nrow(anova_df) > 0) {
    
    # Check if the CSV file already exists
    if (file.exists(anova_csv_path)) {
      existing_anova_results <- tryCatch({
        read.csv(anova_csv_path, stringsAsFactors = FALSE)
      }, error = function(e) {
        # Handle errors during reading
        message("Error reading the existing ANOVA results file: ", e)
        data.frame()
      })
      
      # Merge new results with existing results, updating duplicates
      merged_anova_results <- existing_anova_results %>%
        filter(!(variable %in% anova_df$variable & term %in% anova_df$term)) %>%
        bind_rows(anova_df)
    } else {
      # If no existing file, use the new results
      merged_anova_results <- anova_df
    }
    
    # Write the merged ANOVA results to CSV and Excel files
    write.csv(merged_anova_results, anova_csv_path, row.names = FALSE)
    write.xlsx(merged_anova_results, anova_xlsx_path)
  } else {
    message("No valid ANOVA test results were generated.")
  }
  
  # Save the Tukey's HSD results
  if (nrow(tukey_df) > 0) {
    
    # Check if the CSV file already exists
    if (file.exists(tukey_csv_path)) {
      existing_tukey_results <- tryCatch({
        read.csv(tukey_csv_path, stringsAsFactors = FALSE)
      }, error = function(e) {
        # Handle errors during reading
        message("Error reading the existing Tukey's HSD results file: ", e)
        data.frame()
      })
      
      # Merge new results with existing results, updating duplicates
      merged_tukey_results <- existing_tukey_results %>%
        filter(!(variable %in% tukey_df$variable & comparison %in% tukey_df$comparison)) %>%
        bind_rows(tukey_df)
    } else {
      # If no existing file, use the new results
      merged_tukey_results <- tukey_df
    }
    
    # Write the merged Tukey's HSD results to CSV and Excel files
    write.csv(merged_tukey_results, tukey_csv_path, row.names = FALSE)
    write.xlsx(merged_tukey_results, tukey_xlsx_path)
  } else {
    message("No valid Tukey's HSD test results were generated.")
  }
  
  # Return the ANOVA and Tukey's HSD result objects for further analysis
  return(list(anova_result = anova_result, tukey_result = tukey_result))
}



