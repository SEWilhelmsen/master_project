# Analyse fibrosis data 
# silje Wilhelmsen


# Preparation if necessary
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/load_libraries_hypothesis_testing.R")
source("C:/Users/siljeew/Master_project/Phenotypic_data/Analysis/prepare_fibrosis_data.R") # Output: fibrosis_data
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_kruskal_wallis_test.R")

# Load libraries
load_libraries()

# Load data and set parameters
data <- fibrosis_data
output_dir_plot <- "C:/Users/siljeew/Master_project/Phenotypic_data/Plots/"
output_dir <- "C:/Users/siljeew/Master_project/Phenotypic_data/Data/"

variable_of_interest <- 'fibrosis'
data$fibrosis <- as.numeric(data$fibrosis)


# Investigate normal distribution
##########################################################################
histogram <- ggplot(data, aes(x = !!sym(variable_of_interest), fill = Condition)) +
  geom_histogram(binwidth = 0.3, color = "lightgrey", alpha = 0.7, position = "identity") +
  labs(title = paste("Histogram of", variable_of_interest),
       x = variable_of_interest,
       y = "Frequency") +
  facet_wrap(Condition ~ Timepoint, ncol = 6) +
  scale_fill_manual(values = c("SHAM" = "grey22", "ORAB" = "coral")) +
  theme_minimal()

# Display the histogram
print(histogram)



# Q-q-plot
qq_plot <- ggplot(data, aes(sample = !!sym(variable_of_interest), color = Condition)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = paste("QQ Plot of", variable_of_interest),
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  facet_wrap(~ Timepoint, ncol = 6) +
  scale_color_manual(values = c("SHAM" = "grey22", "ORAB" = "coral")) +
  theme_minimal()

print(qq_plot)


# Get n in each group
##########################################################
sample_counts <- data %>%
  group_by(Timepoint, Condition) %>%
  summarize(sample_count = n())
print(sample_counts)



# # T-test and wilcoxon test
# ######################################################################
# 
# # Initialize lists to store results
# results_ttest <- list()
# results_wilcox <- list()
# 
# # List of unique time points
# time_points <- unique(data$Timepoint)
# 
# # Loop through each time point for comparisons
# for (tp in time_points) {
#   # Subset data by time point
#   subset_data <- data %>% filter(Timepoint == tp)
#   
#   # Further subset data by condition
#   sham_data <- subset_data %>% filter(Condition == "SHAM") %>% select(variable_of_interest) %>% unlist() # is x
#   orab_data <- subset_data %>% filter(Condition == "ORAB") %>% select(variable_of_interest) %>% unlist() # is y
#   
#   # Perform Student's t-test with explicit equal variance assumption
#   t_test_result <- t.test(sham_data, orab_data, var.equal = TRUE)
#   results_ttest[[tp]] <- t_test_result
#   
#   # Perform Wilcoxon rank-sum test (if non-parametric test is preferred)
#   wilcox_test_result <- wilcox.test(sham_data, orab_data)
#   results_wilcox[[tp]] <- wilcox_test_result
# }
# 
# # Print results
# print("T-test results (Student's):")
# print(results_ttest)
# 
# print("Wilcoxon test results:")
# print(results_wilcox)



# Non-normal data: Kruskal Wallis + Dunn test
############################################################################
output_csv_path <- file.path(output_dir, "kruskal_dunn_fiborisis.csv")
output_xlsx_path <- file.path(output_dir, "kruskal_dunn_fibrosis.xlsx")
processed_data <- data
kruskal_dunn_test_results <- perform_kruskal_dunn(processed_data, variable_of_interest, output_csv_path, output_xlsx_path)

View(kruskal_dunn_test_results)

value <- kruskal_dunn_test_results %>%
  filter(variable == variable_of_interest & group1_timepoint == group2_timepoint) %>%
  select(variable, group1, group2, group1_timepoint, group2_timepoint, p_value_adjusted, median_group1, median_group2 ,iqr_group1, iqr_group2)

print(value)

# Optionally, format output for better presentation in a table
formatted_results <- value %>%
  mutate(
    group_comparison = paste(group1, "vs", group2),
    median_comparison = paste(median_group1, "/", median_group2)
  ) %>%
  select(variable, group1_timepoint, group_comparison, p_value_adjusted, median_comparison)

print(formatted_results)




# ANOVA + Tukey HSD
# table(data$variable_of_interest)
# table(processed_data$variable_of_interest)
# 
# source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_anova_tukey.R")
# processed_data <- data
# colnames(processed_data)
# View(processed_data)
# 
# processed_data$fibrosis <- as.numeric(processed_data$fibrosis)
# processed_data$Condition <- as.factor(processed_data$Condition)
# processed_data$Timepoint <- as.factor(processed_data$Timepoint)
# 
# output_dir <- "C:/Users/siljeew/Master_project/Phenotypic_data/Data"
# 
# anova_csv_path <- file.path(output_dir, "anova_results_fibrosis.csv")
# anova_xlsx_path <- file.path(output_dir, "anova_results_fibrosis.xlsx")
# tukey_csv_path <- file.path(output_dir, "tukey_results_fibrosis.csv")
# tukey_xlsx_path <- file.path(output_dir, "tukey_results_fibrosis.xlsx")
# 
# perform_anova_tukey_test(processed_data, variable_of_interest, anova_csv_path, anova_xlsx_path, tukey_csv_path, tukey_xlsx_path)
# 
file.edit("C:/Users/siljeew/Master_project/Phenotypic_data/Analysis/run_create_boxplot_fibrosis.R")