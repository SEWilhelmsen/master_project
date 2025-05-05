# t-test for progeny scores for one specific pathway by Condition
# Silje Wilhelmsen

perform_pairwise_t_test <- function(data, group_col, activity_col, time_col, pathways_of_interest) {
  t_test_results <- list()
  
  for (pathway in pathways_of_interest) {
    # Subset data for the specific pathway of interest
    pathway_data <- data[data$Pathway == pathway, ]
    
    # Split data by the given timepoint
    split_data <- split(pathway_data, pathway_data[[time_col]])
    
    for (time_point in names(split_data)) {
      data_subset <- split_data[[time_point]]
      
      # Ensure there are enough entries for both groups
      if (nrow(data_subset) > 0) {
        sham_values <- data_subset[[activity_col]][data_subset[[group_col]] == "SHAM"]
        orab_values <- data_subset[[activity_col]][data_subset[[group_col]] == "ORAB"]
        
        if (length(sham_values) >= 2 && length(orab_values) >= 2) {
          # Perform pairwise t-test
          t_test_result <- t.test(orab_values, sham_values, var.equal = TRUE)  # Assuming equal variance; adjust as needed
          
          # Store results
          t_test_results[[length(t_test_results) + 1]] <- data.frame(
            Pathway = pathway,
            Timepoint = time_point,
            p_value = t_test_result$p.value,
            t_statistic = t_test_result$statistic,
            mean_sham = mean(sham_values, na.rm = TRUE),
            mean_orab = mean(orab_values, na.rm = TRUE),
            conf_low = t_test_result$conf.int[1],
            conf_high = t_test_result$conf.int[2],
            stringsAsFactors = FALSE
          )
        } else {
          message(paste("Not enough data for pathway:", pathway, "at Timepoint:", time_point))
        }
      }
    }
  }
  
  # Combine all results into a single data frame
  t_test_results_df <- bind_rows(t_test_results)
  
  # Perform Benjamini-Hochberg correction for p-values
  if (nrow(t_test_results_df) > 0) {
    t_test_results_df$p_adjusted_bh <- p.adjust(t_test_results_df$p_value, method = "BH")
  }
  
  return(t_test_results_df)
}





