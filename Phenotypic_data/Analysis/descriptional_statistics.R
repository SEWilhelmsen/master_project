# Descriptional statistics 
# silje Wilhelmsen


# Source scripts
source("C:/Users/siljeew/Master_project/Phenotypic_data/prepare_fibrosis_data.R") # Output: fibrosis_data
#source("C:/Users/siljeew/Master_project/Phenotypic_data/prepare_cell_size_data.R") # output: cell_size_data_filtered

# Load libraries
library(openxlsx)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyverse)
library(writexl)
library(tidyr)
library(car)
library(summarytools)

# Load data and set parameters
# If not created cell size data, load the saved file
#data <- read_excel("Data/cell_size_data_filtered.xlsx")
data <- cell_size_data_filtered

# data <- fibrosis_data
View(data)

output_dir_plot <- "C:/Users/siljeew/Master_project/Phenotypic_data/Plots/"

variable_of_interest <- 'width_um'
data$`width_um` <- as.numeric(data$`width_um`)
#data <- data %>% filter(area_um2 < 10000)
#data <- data %>% filter(length_um < 1000)
   
data <- data %>% filter(width_um < 100)


# Investigate normal distribution
##########################################################################
histogram <- ggplot(data, aes(x = !!sym(variable_of_interest), fill = Condition)) +
  geom_histogram(binwidth = 0.5, color = "black", alpha = 0.7, position = "identity") + #For cell size set binwidth to 500
  labs(title = paste("Histogram of", variable_of_interest),
       x = variable_of_interest,
       y = "Frequency") +
  facet_wrap(~ Timepoint, ncol = 6) +
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

