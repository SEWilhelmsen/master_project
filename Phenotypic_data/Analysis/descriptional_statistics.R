# Descriptional statistics 
# silje Wilhelmsen



# Preparation
###############################################################################
source("C:/Users/siljeew/Master_project/Phenotypic_data/Analysis/load_libraries_for_plots_phenotype_data.R")
load_libraries()


#source("C:/Users/siljeew/Master_project/Phenotypic_data/Analysis/prepare_fibrosis_data.R") # Output: fibrosis_data
# Preparational scripts
#source("C:/Users/siljeew/Master_project/Phenotypic_data/Analysis/prepare_cell_size_data.R") # output: cell_size_data_filtered


# Load data and set parameters
# If not created cell size data, load the saved file
cell_size_data_filtered <- read_excel("C:/Users/siljeew/Master_project/Phenotypic_data/Data/cell_size_data_filtered.xlsx")
data <- cell_size_data_filtered

# data <- fibrosis_data
View(data)

output_dir_plot <- "C:/Users/siljeew/Master_project/Phenotypic_data/Plots/"

variable_of_interest <- 'width_um'
data$`width_um` <- as.numeric(data$`width_um`)
data <- data %>% filter(width_um < 100)

# variable_of_interest <- 'area_um2'
# data$`area_um2` <- as.numeric(data$`area_um2`)
# data <- data %>% filter(area_um2 < 10000)


# colnames(data)
# variable_of_interest <- 'length_um'
# data$length_um <- as.numeric(data$length_um)
# data <- data %>% filter(length_um < 1000)



# Investigate normal distribution
##########################################################################
histogram <- ggplot(data, aes(x = !!sym(variable_of_interest), fill = Condition)) +
  geom_histogram(binwidth = 10, color = "lightgrey", alpha = 0.7, position = "identity") + #For cell size set binwidth to 500
  labs(title = paste("Histogram of", variable_of_interest),
       x = variable_of_interest,
       y = "Frequency") +
  facet_wrap(Condition ~ Timepoint, ncol = 6) +
  scale_fill_manual(values = c("SHAM" = "grey22", "ORAB" = "coral")) +
  theme_minimal()

# Display the histogram
print(histogram)

# Width
histogram_width <- ggplot(data, aes(x = !!sym(variable_of_interest), fill = Condition)) +
  geom_histogram(binwidth = 10, color = "lightgrey", alpha = 0.7, position = "identity") + #For cell size set binwidth to 500
  labs(title = paste("Histogram of", variable_of_interest),
       x = variable_of_interest,
       y = "Frequency") +
  facet_wrap(Condition ~ Timepoint, ncol = 6) +
  scale_fill_manual(values = c("SHAM" = "grey22", "ORAB" = "coral")) +
  theme_minimal()

# Display the histogram
print(histogram_width)




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

# Display the QQ plot
print(qq_plot)


# Get n in each group
##########################################################
sample_counts <- data %>%
  group_by(Timepoint, Condition) %>%
  summarize(sample_count = n())

# View the result
print(sample_counts)

# Find number of animals
animal_counts <- data %>%
  group_by(Timepoint, Condition) %>%
  summarize(unique_animals = n_distinct(sample_id))

print(animal_counts)


# T-test and wilcoxon test
######################################################################

# Initialize lists to store results
results_ttest <- list()
results_wilcox <- list()

# List of unique time points
time_points <- unique(data$Timepoint)

# Loop through each time point for comparisons
for (tp in time_points) {
  # Subset data by time point
  subset_data <- data %>% filter(Timepoint == tp)

  # Further subset data by condition
  sham_data <- subset_data %>% filter(Condition == "SHAM") %>% select(variable_of_interest) %>% unlist() # is x
  orab_data <- subset_data %>% filter(Condition == "ORAB") %>% select(variable_of_interest) %>% unlist() # is y

  # Perform Student's t-test with explicit equal variance assumption
  t_test_result <- t.test(sham_data, orab_data, var.equal = TRUE)
  results_ttest[[tp]] <- t_test_result

  # Perform Wilcoxon rank-sum test (if non-parametric test is preferred)
  wilcox_test_result <- wilcox.test(sham_data, orab_data)
  results_wilcox[[tp]] <- wilcox_test_result
}

# Print results
print("T-test results (Student's):")
print(results_ttest)

print("Wilcoxon test results:")
print(results_wilcox)


# ANOVA + Tukey HSD
###########################################################################
table(data$variable_of_interest)

table(processed_data$variable_of_interest)

source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_anova_tukey.R")
processed_data <- data
colnames(processed_data)
View(processed_data)

processed_data$length_um <- as.numeric(processed_data$length_um)
processed_data$Condition <- as.factor(processed_data$Condition)
processed_data$Timepoint <- as.factor(processed_data$Timepoint)

output_dir <- "C:/Users/siljeew/Master_project/Phenotypic_data/Data"

anova_csv_path <- file.path(output_dir, "anova_results_cell_size.csv")
anova_xlsx_path <- file.path(output_dir, "anova_results_cell_size.xlsx")
tukey_csv_path <- file.path(output_dir, "tukey_results_cell_size.csv")
tukey_xlsx_path <- file.path(output_dir, "tukey_results_cell_size.xlsx")

perform_anova_tukey_test(processed_data, variable_of_interest, anova_csv_path, anova_xlsx_path, tukey_csv_path, tukey_xlsx_path)



data$Timepoint <- as.factor(data$Timepoint)

ggplot(data, aes(x = Timepoint, y = length_um, fill = Condition)) +
  geom_violin()



################################################################################
# Investigate descriptional statistics
freq(data)

# summary(data)
# t.test(data$Condition)

# Calculate summary statistics for each condition

# Calculate the confidence interval
result <- t.test(data$width_um)
# Extract the confidence interval
confidence_interval <- result$conf.int
confidence_interval

