# Function to perform ANOVA and Tukey HSD for one gene by stress status
# Silje Wilhelmsen


# Load required libraries
library(dplyr)
library(readr) 
library(openxlsx) 

# perform_anova_bh_for_gene <- function(data, Stress_Status, expression, Timepoint, gene_of_interest, output_csv_path, output_xlsx_path) {
#   # Ensure column names exist and are strings
#   stopifnot(is.character(Stress_Status), is.character(expression), is.character(Timepoint))
#   stopifnot(all(c(Stress_Status, expression, Timepoint) %in% names(data)))
# 
#   # Subset data for the specific gene of interest
#   gene_data <- data[data$gene == gene_of_interest, ]
# 
#   # Debugging step: Print the gene subset outcome
#   if (nrow(gene_data) == 0) {
#     stop(paste("No data found for gene:", gene_of_interest))
#   } else {
#     print("Subsetting successful:")
#     print(head(gene_data))
#   }
# 
#   # Split data by the given timepoint
#   split_data <- split(gene_data, gene_data[[Timepoint]])
#   results <- list()
#   comparisons_results <- list()  # List to store pairwise comparison results
# 
#   # Conduct ANOVA for each timepoint
#   for (time_point in names(split_data)) {
#     data_subset <- split_data[[time_point]]
#     print(paste("Processing Timepoint:", time_point))
#     print(head(data_subset))  # Added visibility for subset
# 
#     if (nrow(data_subset) >= 3) {  # Ensure enough data for ANOVA
#       # Ensure Stress_Status is a factor
#       data_subset[[Stress_Status]] <- factor(data_subset[[Stress_Status]])
# 
#       # Perform ANOVA
#       anova_result <- aov(as.formula(paste(expression, "~", Stress_Status)), data = data_subset)
#       print(summary(anova_result))  # Debug ANOVA summary
# 
#       # Extract p-value from ANOVA summary
#       p_value <- summary(anova_result)[[1]][["Pr(>F)"]][1]
# 
#       # Store the p-value
#       results[[time_point]] <- data.frame(
#         gene = gene_of_interest,
#         Timepoint = time_point,
#         p_value = p_value
#       )
# 
#       # Perform pairwise comparisons if there are enough unique groups
#       if (length(unique(data_subset[[Stress_Status]])) > 1) {  # Check if more than one group exists
#         pairwise_results <- pairwise.t.test(data_subset[[expression]], 
#                                             data_subset[[Stress_Status]], 
#                                             p.adjust.method = "BH", 
#                                             pool.sd = FALSE)  # Setting pool.sd to FALSE to perform Welch's t-test
#         
#         # Convert pairwise results to a data frame
#         pairwise_df <- as.data.frame(pairwise_results$p.value)
#         
#         if (!is.null(pairwise_df)) {  # Check if the resulting data frame exists
#           pairwise_df <- as.data.frame(as.table(pairwise_df))  # Convert to long format
#           pairwise_df <- pairwise_df[!is.na(pairwise_df$p.value), ]  # Remove NA comparisons
#           
#           # Rename grouping columns
#           colnames(pairwise_df) <- c("group1", "group2", "p_value")
#           
#           # Add gene and timepoint to pairwise results
#           pairwise_df$gene <- gene_of_interest
#           pairwise_df$Timepoint <- time_point
#           
#           # Store pairwise comparison results
#           comparisons_results[[time_point]] <- pairwise_df
#         } else {
#           warning(paste("No comparisons can be made for Timepoint:", time_point))
#         }
#       } else {
#         warning(paste("Not enough unique groups for pairwise comparison at Timepoint:", time_point))
#       }
#     }
#   }
# 
#   results_df <- bind_rows(results)
# 
#   # Apply Benjamini-Hochberg correction for single ANOVA results
#   results_df <- results_df %>%
#     mutate(p_adj = p.adjust(p_value, method = "BH"))
# 
#   comparisons_results_df <- bind_rows(comparisons_results)
# 
#   # Update existing results
#   if (file.exists(output_csv_path)) {
#     existing_results <- read_csv(output_csv_path, show_col_types = FALSE)
# 
#     # Remove old entries for the current gene and timepoints, then merge
#     updated_results <- existing_results %>%
#       filter(!(gene == gene_of_interest & Timepoint %in% results_df$Timepoint))
# 
#     results_df <- bind_rows(updated_results, results_df)
#   }
# 
#   # Write updated results
#   write_csv(results_df, output_csv_path)
#   write.xlsx(results_df, output_xlsx_path, rowNames = FALSE)
# 
#   # Write pairwise comparison results
#   pairwise_output_csv_path <- gsub(".csv", "_pairwise.csv", output_csv_path)
#   pairwise_output_xlsx_path <- gsub(".xlsx", "_pairwise.xlsx", output_xlsx_path)
# 
#   write_csv(comparisons_results_df, pairwise_output_csv_path)
#   write.xlsx(comparisons_results_df, pairwise_output_xlsx_path, rowNames = FALSE)
# 
#   return(list(anova_results = results_df, pairwise_comparisons = comparisons_results_df))
# }
# 
# 
# anova_bh_gene_results <- perform_anova_bh_for_gene(
#   data = combined_data_long,
#   Stress_Status = "Stress_Status",
#   expression = "expression",
#   Timepoint = "Timepoint",
#   gene_of_interest = gene_of_interest,
#   output_csv_path = "C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/anova_bh_gene_results.csv",
#   output_xlsx_path = "C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/anova_bh_gene_results.xlsx"
# )


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



