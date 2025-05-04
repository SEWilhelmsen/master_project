# Hypothesis testing
# Input: dataset of small size (e.g. n=10 per group)
# Output: test for normality, histogram, qq-plot, hypothesis test (Kruskal Wallis, )
# Silje Wilhelmsen

# Preparation
###################################################################################3

# Source necessary files
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/load_libraries_hypothesis_testing.R")
source("C:/Users/siljeew/Master_project/Phenotypic_data/Analysis/extract_and_clean_data_from_animal_overview.R")
source("C:/Users/siljeew/Master_project/Phenotypic_data/Analysis/perform_shapiro_wilks_test.R")
source("C:/Users/siljeew/Master_project/Phenotypic_data/Analysis/create_qqplot_histogram_weight.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_kruskal_wallis_test.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_t_test.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_welch_test.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_anova_tukey.R")




# Load libraries
load_libraries()

# File path to the Excel file
file_path <- "C:/Users/siljeew/Master_project/Phenotypic_data/Data/Animal_overview_TP.xlsx"


# Prepare data from original excel file
processed_data <- process_excel_data(file_path, output_view = TRUE) 
View(processed_data) # Check the results

# Define variables of interest and output directory
variable_of_interest <- 'lvw_bw'  # Change this to the desired variable
output_dir <- "C:/Users/siljeew/Master_project/Phenotypic_data/Data"

# Ensure the directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}



# Assess normality
####################
# Perform Shapiro-Wilk test
output_csv_path <- file.path(output_dir, "shapiro_wilks_results.csv")
perform_shapiro_wilks_test(processed_data, variable_of_interest, output_csv_path)


# Visual investigation: qq-plot
output_dir_distribution <- "C:/Users/siljeew/Master_project/Phenotypic_data/Plots/Distribution"
qq_plots <- create_qqplot(processed_data, variable_of_interest, output_dir_distribution)

# Visual investigation: histogram
# You can change the binwidth to adjust for level of detail. 
histograms <- create_histogram(processed_data, variable_of_interest, output_dir_distribution)

histogram_plot <- ggplot(processed_data, aes(x = 'lw', fill = Condition, stat = count)) +
  geom_histogram(binwidth = 2, color = "black", alpha = 0.7, position = "dodge") +
  facet_grid(Condition ~ Timepoint) +
  labs(title = "", 
       x = 'LW (mg)', 
       y = "Frequency") +
  scale_fill_manual(values = c("sham" = "grey22", "ORAB" = "coral")) + 
  theme_minimal() +
  theme(legend.position = "none")

# Print the plot
print(histogram_plot)



colnames(processed_data)

lung_data <- processed_data %>%
  rename(lw = "LW_(mg)")
histogram_lw <- lung_data %>%
  ggplot(aes(x = lw, fill = Condition)) +
  geom_histogram(binwidth= 3, color = "#e9ecef", alpha=1, position = 'identity') +
  scale_fill_manual(values = c("SHAM" = "grey22", "ORAB" = "coral")) +
  theme_classic() +
  facet_grid(Condition ~ Timepoint)

histogram_lw

processed_data <- processed_data %>%
  rename(lw = "LW_(mg)")
histogram_hw <- processed_data %>%
  ggplot(aes(x = hw, fill = Condition)) +
  geom_histogram(binwidth= 6, color = "#e9ecef", alpha=1, position = 'identity') +
  scale_fill_manual(values = c("SHAM" = "grey22", "ORAB" = "coral")) +
  theme_classic() +
  facet_grid(Condition ~ Timepoint)

histogram_hw

histogram_bw <- processed_data %>%
  ggplot(aes(x = bw, fill = Condition)) +
  geom_histogram(binwidth= 0.7, color = "#e9ecef", alpha=1, position = 'identity') +
  scale_fill_manual(values = c("SHAM" = "grey22", "ORAB" = "coral")) +
  theme_classic() +
  facet_grid(Condition ~ Timepoint)

histogram_bw


# Hypothesis testing
######################
# Kruskal Wallis test for significance between multiple groups 
output_csv_path <- file.path(output_dir, "kruskal_dunn_results.csv")
output_xlsx_path <- file.path(output_dir, "kruskal_dunn_results.xlsx")
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



# Mann-Whitney U test for significance between two independent groups,
output_csv_path <- file.path(output_dir, "mann_whitney_results.csv")
output_xlsx_path <- file.path(output_dir, "mann_whitney_results.xlsx")
mann_whitney_results <- perform_mann_whitney_test(processed_data, variable_of_interest, output_csv_path, output_xlsx_path)


# t-test for parametric samples
output_csv_path <- file.path(output_dir, "t_test_results.csv")
output_xlsx_path <- file.path(output_dir, "t_test_results.xlsx")
perform_t_test(processed_data, variable_of_interest, output_csv_path, output_xlsx_path)


# Perform Welch test between two groups with unequal variance
output_csv_path <- file.path(output_dir, "welch_test_results.csv")
output_xlsx_path <- file.path(output_dir, "welch_test_results.xlsx")
perform_welch_test(processed_data, variable_of_interest, output_csv_path, output_xlsx_path)


# ANOVA and TukeyHSD
processed_data <- processed_data %>%
  rename(lw = "LW_(mg)")
variable_of_interest <- 'lw'
View(processed_data)

anova_csv_path <- file.path(output_dir, "anova_results.csv")
anova_xlsx_path <- file.path(output_dir, "anova_results.xlsx")
tukey_csv_path <- file.path(output_dir, "tukey_results.csv")
tukey_xlsx_path <- file.path(output_dir, "tukey_results.xlsx")

perform_anova_tukey_test(processed_data, variable_of_interest, anova_csv_path, anova_xlsx_path, tukey_csv_path, tukey_xlsx_path)


# Combine Condition and Timepoint into a group identifier
processed_data <- processed_data %>%
  mutate(Group_Timepoint = paste(Condition, Timepoint, sep = "_"))

# Correction method: none
pairwise_results <- pairwise.t.test(processed_data$lw, processed_data$Group_Timepoint, p.adjust.method = "none", alternative = c("two.sided"))
print(pairwise_results)

# Correction method: Benjamini-Hochberg
pairwise_results <- pairwise.t.test(processed_data$lw, processed_data$Group_Timepoint, p.adjust.method = "BH", alternative = c("two.sided"))
print(pairwise_results)


# Visualization
# file.edit("C:/Users/siljeew/Master_project/Phenotypic_data/Analysis/create_boxplot_weight.R")


