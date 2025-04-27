# Kruskal wallis test 
# Non parametric test of multiple groups. 
# H0: the groups are equal 
# If significant (AKA groups are not equal), perform Dunn test to test all possible pairs
# or perform with no correction method to investigate crude p-value 
# or use method = "bh" for Benjamini Hochberg correction


library(dplyr)
library(dunn.test)
library(openxlsx) 

perform_kruskal_dunn <- function(data, variable_of_interest, output_csv_path, output_xlsx_path) {
  # Combine Condition and Timepoint into a group identifier
  data <- data %>%
    mutate(group = paste(Condition, Timepoint, sep = "_"))
  
<<<<<<< HEAD
  # Perform Kruskal-Wallis test
  kruskal_result <- kruskal.test(as.formula(paste(variable_of_interest, "~ group")), data = data)
=======
  # Create a new column 'group' by combining 'condition' and 'time_point'
  processed_data <- processed_data %>%
    mutate(group = paste(Condition, Timepoint, sep = "_"))
>>>>>>> 7df848ad59037d1ccdbeb933fda380e63ccfea21
  
  if (kruskal_result$p.value < 0.05) {
    # If Kruskal-Wallis is significant, perform post-hoc Dunn's test
    # dunn_results <- dunn.test(data[[variable_of_interest]], data$group, method = "bonferroni")
    # dunn_results <- dunn.test(data[[variable_of_interest]], data$group, method = "none")
    dunn_results <- dunn.test(data[[variable_of_interest]], data$group, method = "bh")
    
    # Extract Dunn's test results into a data frame
    result_pairs <- dunn_results$comparisons
    p_values_adjusted <- dunn_results$P.adjusted
    
    # Calculate medians and IQRs for each comparison group
    group_stats <- data %>%
      group_by(group) %>%
      summarize(
        median = median(get(variable_of_interest), na.rm = TRUE),
        iqr = IQR(get(variable_of_interest), na.rm = TRUE)
      )
    
    # Split comparison pairs into group1 and group2
    comparison_split <- strsplit(result_pairs, " - ")
    comparison_df <- do.call(rbind, comparison_split)
    colnames(comparison_df) <- c("group1", "group2")
    
    # Extract timepoints from group identifiers
    group1_timepoints <- sub(".*_", "", comparison_df[, "group1"])
    group2_timepoints <- sub(".*_", "", comparison_df[, "group2"])
    
    # Merge median and IQR with comparisons
    stats_group1 <- group_stats[group_stats$group %in% comparison_df[, "group1"], ]
    stats_group2 <- group_stats[group_stats$group %in% comparison_df[, "group2"], ]
    
    # Create result data frame with median and IQR
    result_df <- data.frame(
      variable = variable_of_interest,
      p_value_adjusted = p_values_adjusted,
      group1 = comparison_df[, "group1"],
      group1_timepoint = group1_timepoints,
      group2 = comparison_df[, "group2"],
      group2_timepoint = group2_timepoints,
      median_group1 = stats_group1$median,
      iqr_group1 = stats_group1$iqr,
      median_group2 = stats_group2$median,
      iqr_group2 = stats_group2$iqr,
      stringsAsFactors = FALSE
    )
    
  } else {
    result_df <- data.frame(
      variable = variable_of_interest,
      p_value_adjusted = NA,
      group1 = NA,
      group1_timepoint = NA,
      group2 = NA,
      group2_timepoint = NA,
      median_group1 = NA,
      iqr_group1 = NA,
      median_group2 = NA,
      iqr_group2 = NA,
      stringsAsFactors = FALSE
    )
  }
  
  # Update existing results
  if (file.exists(output_csv_path)) {
    existing_results <- read.csv(output_csv_path, stringsAsFactors = FALSE)
    # Remove old entries for the current variable_of_interest, then merge
    updated_results <- existing_results %>%
      filter(!(variable == variable_of_interest))
    
    result_df <- bind_rows(updated_results, result_df)
  }
  
  # Write updated result to file
  write.csv(result_df, output_csv_path, row.names = FALSE)
  write.xlsx(result_df, output_xlsx_path)
  
  return(result_df)
}


 
# # Function to perform Kruskal-Wallis test and save the results
# perform_kruskal_wallis_test <- function(processed_data, variable, output_path) {
#   
#   # Create a new column 'group' by combining 'condition' and 'time_point'
#   processed_data <- processed_data %>%
#     mutate(group = paste(Condition, Timepoint, sep = "_"))
#   
#   result <- kruskal.test(as.formula(paste(variable, "~ group")), data = processed_data)
#   
#   # Convert the htest object to a data frame
#   result_df <- data.frame(
#     variable_of_interest = variable,
#     statistic = result$statistic,
#     parameter = result$parameter,
#     p_value = result$p.value,
#     method = result$method,
#     data_name = result$data.name,
#     stringsAsFactors = FALSE
#   )
#   
#   View(result_df)
#   
#   # Check if the CSV file already exists
#   if (file.exists(output_path)) {
#     existing_results <- read.csv(output_path)
#     
#     # Merge new results with existing results, updating duplicates
#     merged_results <- existing_results %>%
#       filter(!(variable_of_interest %in% result_df$variable_of_interest)) %>%
#       bind_rows(result_df)
#   } else {
#     # If no existing file, use the new results
#     merged_results <- result_df
#   }
#   
#   # Save the results to a CSV file
#   write.csv(merged_results, output_path, row.names = FALSE)
#   
#   View(merged_results)
#   
#   # Optional: View the results
#   return(merged_results)
# }


