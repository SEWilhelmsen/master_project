# t-test for progeny scores for one specific pathway
# Silje Wilhelmsen

group1 <- "SHAM CM"
# group2 <- "Not stressed CM"
group2 <- "Stressed CM"

# Prepare a list to store the results
t_test_results <- list()

# Perform t-tests for each pathway at each timepoint
for (pathway in unique(progeny_results_filtered$Pathway)) {
  # Subset data for the current pathway
  pathway_data <- progeny_results_filtered %>%
    filter(Pathway == pathway, Stress_Status %in% c(group1, group2))
  
  # Split by Timepoint and perform t-tests
  for (timepoint in unique(pathway_data$Timepoint)) {
    data_subset <- pathway_data %>%
      filter(Timepoint == timepoint)
    
    # Ensure there are enough observations per group
    group1_values <- data_subset %>% filter(Stress_Status == group1) %>% pull(Activity)
    group2_values <- data_subset %>% filter(Stress_Status == group2) %>% pull(Activity)
    
    if (length(group1_values) >= 2 && length(group2_values) >= 2) {
      # Perform the t-test
      t_test_result <- t.test(group2_values, group1_values)
      
      # Store results
      t_test_results[[length(t_test_results) + 1]] <- data.frame(
        Pathway = pathway,
        Timepoint = timepoint,
        p_value = t_test_result$p.value,
        t_statistic = t_test_result$statistic,
        group2 = group2,
        mean_group2 = mean(group2_values, na.rm = TRUE),
        group1 = group1,
        mean_group1 = mean(group1_values, na.rm = TRUE),
        conf_low = t_test_result$conf.int[1],
        conf_high = t_test_result$conf.int[2],
        stringsAsFactors = FALSE
      )
    } else {
      message(paste("Not enough data for pathway:", pathway, "at timepoint:", timepoint))
    }
  }
}

# Combine all results into a single data frame
t_test_results_df <- bind_rows(t_test_results)

# Perform Benjamini-Hochberg correction for p-values
if (nrow(t_test_results_df) > 0) {
  t_test_results_df$p_adjusted_bh <- p.adjust(t_test_results_df$p_value, method = "BH")
}

View(t_test_results_df)
# View results
print(t_test_results_df)


# Add label
t_test_results_df <- t_test_results_df %>%
  mutate(p_adjusted_label_bh = ifelse(p_adjusted_bh < 0.001, "***",
                                      ifelse(p_adjusted_bh < 0.01, "**",
                                             ifelse(p_adjusted_bh < 0.05, "*",
                                                    "ns"))))


