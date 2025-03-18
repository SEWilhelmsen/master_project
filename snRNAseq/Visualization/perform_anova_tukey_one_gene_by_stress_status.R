# Function to perform ANOVA and Tukey HSD for one gene by stress status
# Silje Wilhelmsen

# Function to perform ANOVA and Tukey HSD per gene
perform_anova_tukey_for_gene <- function(data, Stress_Status, expression, Timepoint, gene_of_interest, output_csv_path, output_xlsx_path) {
  # Ensure column names exist and are strings
  stopifnot(is.character(Stress_Status), is.character(expression), is.character(Timepoint))
  stopifnot(all(c(Stress_Status, expression, Timepoint) %in% names(data)))
  
  # Subset data for the specific gene of interest
  gene_data <- data[data$gene == gene_of_interest, ]
  
  # Debugging step: Print the gene subset outcome
  if (nrow(gene_data) == 0) {
    stop(paste("No data found for gene:", gene_of_interest))
  } else {
    print("Subsetting successful:")
    print(head(gene_data))
  }
  
  # Split data by the given timepoint
  split_data <- split(gene_data, gene_data[[Timepoint]])
  results <- list()
  
  # Conduct ANOVA and Tukey's test
  for (time_point in names(split_data)) {
    data_subset <- split_data[[time_point]]
    print(paste("Processing Timepoint:", time_point))
    print(head(data_subset))  # Added visibility for subset
    
    if (nrow(data_subset) >= 3) {  # Ensure enough data for ANOVA
      # Ensure Stress_Status is a factor
      data_subset[[Stress_Status]] <- factor(data_subset[[Stress_Status]])
      
      # Perform ANOVA
      anova_result <- aov(expression ~ Stress_Status, data = data_subset)
      print(summary(anova_result))  # Debug ANOVA summary
      
      # Tukey HSD post-hoc testing
      tukey_result <- TukeyHSD(anova_result)
      
      # Access the correct matrix from Tukey result
      tukey_df <- as.data.frame(tukey_result[["Stress_Status"]])
      
      # Add rownames as a column for group comparisons
      tukey_df$comparison <- rownames(tukey_df)
      
      # Add group1 and group2 derived from split parts of comparison
      comparison_parts <- strsplit(tukey_df$comparison, "-")
      tukey_df$group1 <- sapply(comparison_parts, `[`, 1)
      tukey_df$group2 <- sapply(comparison_parts, `[`, 2)
      
      # Add gene and timepoint to results
      tukey_df$gene <- gene_of_interest
      tukey_df$Timepoint <- time_point
      
      # Collect results
      results[[time_point]] <- tukey_df
    }
  }
  
  results_df <- bind_rows(results)
  
  # Update existing results
  if (file.exists(output_csv_path)) {
    existing_results <- read_csv(output_csv_path, show_col_types = FALSE)
    
    # Remove old entries for the current gene and timepoints, then merge
    updated_results <- existing_results %>%
      filter(!(gene == gene_of_interest & Timepoint %in% results_df$Timepoint))
    
    results_df <- bind_rows(updated_results, results_df)
  }
  
  # Write updated results
  write_csv(results_df, output_csv_path)
  write.xlsx(results_df, output_xlsx_path)
  
  return(results_df)
}




