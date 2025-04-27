# perform t-test on a gene ontology process or another sum of genes
# Silje Wilhelmsen

# Load libraries
library(dplyr)
library(readr)
library(writexl)


# Function to perform t-tests and update existing results
perform_t_test <- function(data, group_col, value_col, time_col, go_process_of_interest, output_csv_path, output_xlsx_path) {
  split_data <- split(data, data[[time_col]])
  results <- list()
  
  # Conduct t-tests
  for (time_point in names(split_data)) {
    data_subset <- split_data[[time_point]]
    
    sham_values <- data_subset %>% filter(!!sym(group_col) == "SHAM") %>% pull(!!sym(value_col))
<<<<<<< HEAD
    orab_values <- data_subset %>% filter(!!sym(group_col) == "ORAB") %>% pull(!!sym(value_col))
=======
    ab_values <- data_subset %>% filter(!!sym(group_col) == "ORAB") %>% pull(!!sym(value_col))
>>>>>>> 7df848ad59037d1ccdbeb933fda380e63ccfea21
    
    if (length(sham_values) >= 2 & length(orab_values) >= 2) {
      t_test_result <- t.test(sham_values, orab_values)
      
      results[[time_point]] <- data.frame(
        process = go_process_of_interest,
        time_point = time_point,
        group1 = "SHAM",
        group2 = "ORAB",
        t_statistic = t_test_result$statistic,
        p_value = t_test_result$p.value,
        conf_low = t_test_result$conf.int[1],
        conf_high = t_test_result$conf.int[2],
        mean_sham = mean(sham_values, na.rm = TRUE),
        mean_orab = mean(orab_values, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  
  results_df <- bind_rows(results)
  
  # Perform Benjamini-Hochberg adjustment for p-values
  results_df$p_adj <- p.adjust(results_df$p_value, method = "BH")
  
  # Update existing results
  if (file.exists(output_csv_path)) {
    existing_results <- read_csv(output_csv_path, show_col_types = FALSE)
    
    # Remove old entries for the current process and time points, then merge
    updated_results <- existing_results %>%
      filter(!(process == go_process_of_interest & time_point %in% results_df$time_point))
    
    results_df <- bind_rows(updated_results, results_df)
  }
  
  # Write updated results
  write_csv(results_df, output_csv_path)
  write_xlsx(results_df, output_xlsx_path)
  
  return(results_df)
}
