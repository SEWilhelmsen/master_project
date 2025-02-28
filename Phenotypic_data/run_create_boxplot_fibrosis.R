# Boxplot 
# Input: dataset of small size (e.g. n=10 per group)
# Output: 
# Silje Wilhelmsen

# Preparation
###############################################################################


# Load libraries
library(openxlsx)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyverse)
library(writexl)
library(tidyr)

# File path to the Excel file
fibrosis_data <- read.xlsx("C:/Users/siljeew/Master_project/Phenotypic_data/Fibrosis.xlsx")
output_dir_plot <- "C:/Users/siljeew/Master_project/Phenotypic_data/Plots/"

unique(fibrosis_data$timepoint)
fibrosis_data$timepoint <- trimws(fibrosis_data$timepoint)


# Convert to factors and specify the desired order for timepoint
fibrosis_data$timepoint <- factor(fibrosis_data$timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
fibrosis_data$condition <- as.factor(fibrosis_data$condition)

# Check levels of 'timepoint'
levels(fibrosis_data$timepoint)


# Boxplot
##########################################################################################
fibrosis_boxplot <- ggplot(fibrosis_data, aes(x = timepoint, y = fibrosis)) +
  geom_boxplot(aes(fill = condition), position = position_dodge(width = 0.5), alpha = 0.7) +
  labs(title = paste("Fibrosis"),
       x = NULL,
       y = "Fibrosis area of section (%)") +
  scale_fill_manual(name = "Condition", values = c("SHAM" = "white", "ORAB" = "coral"), 
                    labels = c("SHAM" = "SHAM", "ORAB" = "ORAB")) +
  theme_minimal() 

# Display the plot
fibrosis_boxplot


ggsave(file.path(output_dir_plot, paste("fibrosis_boxplot.pdf", sep = "")), plot = fibrosis_boxplot, width = 8, height = 6)


# Barplot
##########################################################################################
ggplot(agg_data, aes(x = timepoint, y = mean_fibrosis, fill = condition)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_fibrosis - se, ymax = mean_fibrosis + se), 
                width = 0.2, position = position_dodge(0.5)) +
  labs(x = NULL, y = "Mean Fibrosis", title = NULL) +
  scale_fill_manual(values = c("SHAM" = "grey", "ORAB" = "coral")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  theme_minimal()


#ggsave(file.path(output_dir_plot, paste("fibrosis_barplot.png", sep = "")), plot = fibrosis_barplot, width = 8, height = 6)
ggsave(file.path(output_dir_plot, paste("fibrosis_barplot.pdf", sep = "")), plot = fibrosis_barplot, width = 8, height = 6)
