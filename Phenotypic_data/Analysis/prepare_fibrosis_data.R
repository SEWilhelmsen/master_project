# Prepare fibrosis data 
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
fibrosis_data <- read.xlsx("C:/Users/siljeew/Master_project/Phenotypic_data/Data/Fibrosis.xlsx")

fibrosis_data <- rename(fibrosis_data, Condition = 'condition', Timepoint = 'timepoint')

unique(fibrosis_data$Timepoint)
fibrosis_data$Timepoint <- trimws(fibrosis_data$Timepoint)

# Convert to factors and specify the desired order for Timepoint
fibrosis_data$Timepoint <- factor(fibrosis_data$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
fibrosis_data$condition <- as.factor(fibrosis_data$condition)

# Check levels of 'Timepoint'
levels(fibrosis_data$Timepoint)

# Output: fibrosis_data