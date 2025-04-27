# Prepare cell size data 
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
cell_size_data <- read_excel("C:/Users/siljeew/Master_project/Phenotypic_data/Cell size statistics.xlsx")
View(cell_size_data)
output_dir_plot <- "C:/Users/siljeew/Master_project/Phenotypic_data/Plots"



# Data preparation
###########################################################################
print(unique(cell_size_data$ID))
colnames(cell_size_data)

# Retrieve specific columns
colnames(cell_size_data)

cell_size_data_filtered <- select(cell_size_data, 'ID', 'Cell-ID', 'Area (um2)', 'Length (um)', 'Width (um)', 'Perimeter (um)')
cell_size_data_filtered <- rename(cell_size_data_filtered, sample_id = 'ID', area_um2 = 'Area (um2)', length_um = 'Length (um)', width_um = 'Width (um)', perimeter_um = 'Perimeter (um)')
View(cell_size_data_filtered)


# Add column of condition and timepoint
# Ikke funnet: "35.5"
sham_samples <- c('60.1', '60.2', '63.1', '63.3', '63.2',
                  '45.1', '47.1', '56.2', '56.1',
                  '41.4', '40.1', '41.1', '41.2', '40.4',
                  '61.1', '61.2', '62.2', '57.1', '57.2',
                  '35.1', '35.3',  '54.3', '55.1', '55.2',
                  '48.2', '59.2', '59.3', '59.4', '59.5')

orab_samples <- c('60.4', '64.1', '64.2', '63.5', '60.3',
                  '45.4', '45.3', '45.5', '56.3', '56.4',
                  '40.2', '40.3', '40.5', '58.1', '58.2',
                  '55.4', '57.4', '57.5', '61.3', '61.4', '62.1',
                  '35.2', '35.4', '55.5', '55.3',
                  '48.3', '48.4', '48.5', '50.1', '50.2')

# Define the conditions for each timepoint
timepoint_06_hours <- c('60.1', '60.2', '63.1', '63.2', '63.3', '60.4', '60.3', '63.5', '64.1', '64.2')
timepoint_12_hours <- c('45.1', '47.1', '56.1', '56.3', '56.2', '56.4', '45.4', '45.3', '45.5')
timepoint_1_day <- c('40.1', '41.4', '40.2', '40.3', '40.4', '40.5', '41.1', '41.2', '58.1', '58.2')
timepoint_3_days <- c('55.4', '57.1', '57.2', '57.4', '57.5', '61.1', '61.2', '61.3', '61.4', '62.1')
timepoint_1_week <- c('35.1', '35.3', '35.2', '35.4', '55.5', '54.3', '55.1', '55.2', '55.3')
timepoint_3_weeks <- c('59.2', '59.3', '59.4', '59.5', '48.2', '48.3', '48.4', '48.5', '50.1', '50.2')

# Create the new 'Timepoint' column and assign values based on the conditions
cell_size_data_filtered$Timepoint <- NA # Initialize column
cell_size_data_filtered$Timepoint <- ifelse(cell_size_data_filtered$sample_id %in% timepoint_06_hours, '6 Hours',
                                            ifelse(cell_size_data_filtered$sample_id %in% timepoint_12_hours, '12 Hours',
                                                   ifelse(cell_size_data_filtered$sample_id %in% timepoint_1_day, '1 Day',
                                                          ifelse(cell_size_data_filtered$sample_id %in% timepoint_3_days, '3 Days',
                                                                 ifelse(cell_size_data_filtered$sample_id %in% timepoint_1_week, '1 Week',
                                                                        '3 Weeks')))))

print(unique(cell_size_data_filtered$Timepoint))

# Create a new column in meta data called Group
cell_size_data_filtered$Condition <- NA
cell_size_data_filtered$Condition <- ifelse(cell_size_data_filtered$sample_id %in% sham_samples, 'SHAM', 'ORAB')
print(unique(cell_size_data_filtered$Condition))


## Create the new group and timepoint
cell_size_data_filtered$Group_Timepoint <- NA
cell_size_data_filtered$Group_Timepoint <- paste(cell_size_data_filtered$Condition, cell_size_data_filtered$Timepoint, sep = "_")
print(unique(cell_size_data_filtered$Group_Timepoint))

# Verify length of data
nrow(cell_size_data)
nrow(cell_size_data_filtered)


# Remove rows with NA in the sample_id column
cell_size_data_filtered <- cell_size_data_filtered[!is.na(cell_size_data_filtered$sample_id), ]
nrow(cell_size_data)
nrow(cell_size_data_filtered)

# Remove rows with NA in the area_um2 column
cell_size_data_filtered <- cell_size_data_filtered[!is.na(cell_size_data_filtered$'area_um2'), ]
nrow(cell_size_data)
nrow(cell_size_data_filtered)


# Save results
View(cell_size_data_filtered)
write.xlsx(cell_size_data_filtered, "C:/Users/siljeew/Master_project/Phenotypic_data/cell_size_data_filtered.xlsx", overwrite=TRUE)


#cell_size_data_filtered <- read.xlsx("C:/Users/Labuser/master_project/Phenotypic_data/cell_size_data_filtered.xlsx")

# Output: cell_size_data_filtered