# Perform ANOVA + TukeyHSD test for gene ontology process
# Silje Wilhelmsen


# Function to perform ANOVA and Tukey's HSD test for a gene ontology process
perform_anova_tukey_for_go_process <- function(data, group_col, value_col, time_col, go_process_of_interest, 
                                               anova_csv_path, anova_xlsx_path, tukey_csv_path, tukey_xlsx_path) {
  # Prepare the formula for ANOVA, including interaction
  formula <- as.formula(paste(value_col, "~", group_col, "*", time_col))
  
  # Perform the ANOVA
  anova_result <- aov(formula, data = data)
  anova_summary <- summary(anova_result)
  
  # Extract ANOVA results into a data frame
  anova_df <- data.frame(
    process = go_process_of_interest,
    term = rownames(anova_summary[[1]]),
    df = anova_summary[[1]]$Df,
    sum_sq = anova_summary[[1]]$"Sum Sq",
    mean_sq = anova_summary[[1]]$"Mean Sq",
    F_value = anova_summary[[1]]$"F value",
    p_value = anova_summary[[1]]$"Pr(>F)",
    stringsAsFactors = FALSE
  )
  
  # Perform Tukey's HSD post-hoc test
  tukey_result <- TukeyHSD(anova_result)
  
  # Extract Tukey's HSD results into a data frame for condition:time_point interaction
  interaction_name <- paste(group_col, time_col, sep = ":")
  
  if (interaction_name %in% names(tukey_result)) {
    tukey_df <- as.data.frame(tukey_result[[interaction_name]], stringsAsFactors = FALSE)
    tukey_df$comparison <- rownames(tukey_df)
    comparison_parts <- strsplit(tukey_df$comparison, "[-:]")
    
    # Create group and timepoint columns
    tukey_df$group1 <- sapply(comparison_parts, `[`, 1)
    tukey_df$group1_timepoint <- sapply(comparison_parts, `[`, 2)
    tukey_df$group2 <- sapply(comparison_parts, `[`, 3)
    tukey_df$group2_timepoint <- sapply(comparison_parts, `[`, 4)
  } else {
    tukey_df <- data.frame()
    message("No interaction term found in TukeyHSD results.")
  }
  
  rownames(tukey_df) <- NULL
  print(tukey_df)
  
  # Rename columns based on observed TukeyHSD structure
  tukey_df <- tukey_df %>%
    rename(
      diff = "diff",
      lwr = "lwr",
      upr = "upr",
      "p adj" = "p adj"  
    ) %>%
    mutate(process = go_process_of_interest, Timepoint = interaction_name) %>%
    dplyr::select(comparison, everything())
  
  # Save the ANOVA results
  update_results(anova_df, anova_csv_path, anova_xlsx_path, c("process", "term"))
  
  # Save the Tukey's HSD results
  update_results(tukey_df, tukey_csv_path, tukey_xlsx_path, c("process", "comparison"))
  
  return(tukey_df)
  
  # Return the ANOVA and Tukey's HSD result data frames
  return(list(anova_df = anova_df, tukey_df = tukey_df))
}

# Helper function to read, update, and write the results
update_results <- function(new_results_df, csv_path, xlsx_path, key_columns) {
  if (nrow(new_results_df) > 0) {
    # Check if the CSV file already exists
    if (file.exists(csv_path)) {
      existing_results <- tryCatch({
        read.csv(csv_path, stringsAsFactors = FALSE)
      }, error = function(e) {
        message("Error reading the existing results file: ", e)
        data.frame()
      })
      
      # Merge new results with existing results, updating duplicates
      merged_results <- existing_results %>%
        filter(!paste0(existing_results[, key_columns], collapse = "-") %in% paste0(new_results_df[, key_columns], collapse = "-")) %>%
        bind_rows(new_results_df)
    } else {
      # If no existing file, use the new results
      merged_results <- new_results_df
    }
    
    # Write the merged results to CSV and Excel files
    write.csv(merged_results, csv_path, row.names = FALSE)
    write.xlsx(merged_results, xlsx_path)
  } else {
    message("No valid test results were generated.")
  }
}


