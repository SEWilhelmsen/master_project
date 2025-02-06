# Hypothesis testing
# Input: dataset of small size (e.g. n=10 per group)
# Output: test for normality, histogram, qq-plot, hypothesis test (Kruskal Wallis, )
# Silje Wilhelmsen

# Preparation
###############

# Source necessary files
source("C:/Users/siljeew/snRNAseq/Analysis/load_libraries_hypothesis_testing.R")
source("C:/Users/siljeew/snRNAseq/Analysis/extract_and_clean_data_from_animal_overview.R")
source("C:/Users/siljeew/snRNAseq/Analysis/perform_shapiro_wilks_test.R")
source("C:/Users/siljeew/snRNAseq/Visualization/create_qqplot_histogram_variable.R")
source("C:/Users/siljeew/snRNAseq/Analysis/perform_kruskal_wallis_test.R")
#source("C:/Users/siljeew/snRNAseq/Analysis/perform_mann_whitney_test.R")
source("C:/Users/siljeew/snRNAseq/Analysis/perform_t_test.R")
source("C:/Users/siljeew/snRNAseq/Analysis/perform_welch_test.R")
source("C:/Users/siljeew/snRNAseq/Analysis/perform_anova_tukey.R")



# Load libraries
load_libraries()

# File path to the Excel file
file_path <- "C:/Users/siljeew/snRNAseq/Animal_overview_TP_copy.xlsx"

# Prepare data from original excel file
processed_data <- process_excel_data(file_path, output_view = TRUE) 
View(processed_data) # Check the results

# Define variables of interest and output directory
variable_of_interest <- 'lvw_bw'  # Change this to the desired variable
output_dir <- "C:/Users/siljeew/snRNAseq/Data"

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
output_dir_distribution <- "Y:/Silje/snRNAseq/Plots/Distribution"
qq_plots <- create_qqplot(processed_data, variable_of_interest, output_dir_distribution)

# Visual investigation: histogram
output_dir_distribution <- "Y:/Silje/snRNAseq/Plots/Distribution"
# You can change the binwidth to adjust for level of detail. 
histograms <- create_histogram(processed_data, variable_of_interest, output_dir_distribution)



# Hypothesis testing
######################

# Kruskal Wallis test for significance between multiple groups 
output_csv_path <- file.path(output_dir, "kruskal_wallis_results.csv")
perform_kruskal_wallis_test(processed_data, variable_of_interest, output_csv_path)


# # Not so relevant with two groups? 
# # Mann-Whitney U test for significance between two independent groups, one variable
# output_csv_path <- file.path(output_dir, "mann_whitney_results.csv")
# perform_mann_whitney_test(processed_data, variable_of_interest, output_csv_path)

# t-test for parametric samples
output_csv_path <- file.path(output_dir, "t_test_results.csv")
output_xlsx_path <- file.path(output_dir, "t_test_results.xlsx")
perform_t_test(processed_data, variable_of_interest, output_csv_path, output_xlsx_path)

# Perform Welch test between two groups with unequal variance
output_csv_path <- file.path(output_dir, "welch_test_results.csv")
output_xlsx_path <- file.path(output_dir, "welch_test_results.xlsx")
perform_welch_test(processed_data, variable_of_interest, output_csv_path, output_xlsx_path)



# ANOVA and TukeyHSD
# ANOVA
anova_csv_path <- file.path(output_dir, "anova_results.csv")
anova_xlsx_path <- file.path(output_dir, "anova_results.xlsx")

# Tukey HSD
tukey_csv_path <- file.path(output_dir, "tukey_results.csv")
tukey_xlsx_path <- file.path(output_dir, "tukey_results.xlsx")

perform_anova_tukey_test(processed_data, variable_of_interest, anova_csv_path, anova_xlsx_path, tukey_csv_path, tukey_xlsx_path)



