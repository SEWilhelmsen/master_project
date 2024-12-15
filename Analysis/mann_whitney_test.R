# Mann-Whitney U test for non-parametric test of difference between two independent groups
# Input:
# Output:
# Silje Wilhelmsen


perform_mann_whitney_test <- function(processed_data, variable_of_interest, output_csv_path) {
  result <- wilcox.test(variable_of_interest, data = processed_data, exact = FALSE)
  
}
