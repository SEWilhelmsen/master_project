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
# library(magrittr)

# File path to the Excel file
cell_size_data <- read_excel("C:/Users/siljeew/Master_project/Cell size statistics.xlsx")
View(cell_size_data)
output_dir_plot <- "C:/Users/siljeew/Master_project/Phenotypic_data/Plots/"

print(unique(cell_size_data$ID))

# Data preparation
###########################################################################
# Use columns "Cell-ID" and "Area (um2)"
colnames(cell_size_data)

cell_size_data_filtered <- select(cell_size_data, 'ID', 'Cell-ID', 'Area (um2)')
cell_size_data_filtered <- rename(cell_size_data_filtered, sample_id = 'ID')




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
timepoint_3_week <- c('59.2', '59.3', '59.4', '59.5', '48.2', '48.3', '48.4', '48.5', '50.1', '50.2')
 
# Create the new 'Timepoint' column and assign values based on the conditions
cell_size_data_filtered$Timepoint <- NA # Initialize column
cell_size_data_filtered$Timepoint <- ifelse(cell_size_data_filtered$sample_id %in% timepoint_1_day, '1 Day',
                                                         ifelse(cell_size_data_filtered$sample_id %in% timepoint_1_week, '1 Week',
                                                                ifelse(cell_size_data_filtered$sample_id %in% timepoint_3_days, '3 Days',
                                                                       ifelse(cell_size_data_filtered$sample_id %in% timepoint_06_hours, '6 Hours',
                                                                              ifelse(cell_size_data_filtered$sample_id %in% timepoint_12_hours, '12 Hours',
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
# Remove rows with NA in the Area (um2) column
cell_size_data_filtered <- cell_size_data_filtered[!is.na(cell_size_data_filtered$'Area (um2)'), ]
nrow(cell_size_data)
nrow(cell_size_data_filtered)

View(cell_size_data_filtered)






# Create graph
#######################################################
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(forcats)

head(cell_size_data_filtered)

cell_size_data_filtered$`Area (um2)` <- as.numeric(cell_size_data_filtered$`Area (um2)`)

cell_size_1w <- cell_size_data_filtered %>%
  filter(Group_Timepoint == "SHAM_1 Week") %>%
  rename(area = `Area (um2)`)

histogram_1w <- ggplot(cell_size_1w, aes(x = area)) +
  geom_histogram(binwidth = 5)

histogram_1w

summary(cell_size_1w)


# Create histogram
histogram <- cell_size_1w %>%
  ggplot(aes(x = area, fill = Group_Timepoint)) +
  geom_histogram(alpha = 0.6, binwidth = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("Area (um²)") +
  ylab("Count") +
  facet_wrap(~ Group_Timepoint)

histogram

cell_size_1w %>%
  ggplot( aes(x=area, fill=Group_Timepoint)) +
  geom_histogram(alpha = 0.6, binwidth = 5) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A histogram wrapping a boxplot") +
  xlab("")

cell_size_1w$area <- as.numeric(cell_size_1w$area)

cell_size_1w %>%
  ggplot(aes(x = area, fill = Group_Timepoint)) +
  geom_histogram(alpha = 0.6, binwidth = 500) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Histogram of Area by Group Timepoint") +
  xlab("Area (um²)")





# Violin plot
# Libraries
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)

# sample size
sample_size = cell_size_1w %>% group_by(Group_Timepoint) %>% summarize(num=n())

# Plot
cell_size_1w %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(Group_Timepoint, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=area, fill=Group_Timepoint)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A Violin wrapping a boxplot") +
  xlab("")


# cell_size_data_filtered$timepoint <- trimws(cell_size_data_filtered$Timepoint)
# 
# 
# # Convert to factors and specify the desired order for timepoint
# cell_size_data$timepoint <- factor(cell_size_data$timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
# cell_size_data$condition <- as.factor(cell_size_data$condition)
# 
# # Check levels of 'timepoint'
# levels(cell_size_data$timepoint)
# 
# 
# # Boxplot
# ##########################################################################################
# fibrosis_boxplot <- ggplot(cell_size_data, aes(x = timepoint, y = fibrosis)) +
#   geom_boxplot(aes(fill = condition), position = position_dodge(width = 0.5), alpha = 0.7) +
#   labs(title = paste("Fibrosis"),
#        x = NULL,
#        y = "Fibrosis area of section (%)") +
#   scale_fill_manual(name = "Condition", values = c("SHAM" = "white", "ORAB" = "coral"), 
#                     labels = c("SHAM" = "SHAM", "ORAB" = "ORAB")) +
#   theme_minimal() 
# 
# # Display the plot
# fibrosis_boxplot
# 
# 
# ggsave(file.path(output_dir_plot, paste("fibrosis_boxplot.pdf", sep = "")), plot = fibrosis_boxplot, width = 8, height = 6)
# 
# 
# # Barplot
# ##########################################################################################
# ggplot(agg_data, aes(x = timepoint, y = mean_fibrosis, fill = condition)) +
#   geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
#   geom_errorbar(aes(ymin = mean_fibrosis - se, ymax = mean_fibrosis + se), 
#                 width = 0.2, position = position_dodge(0.5)) +
#   labs(x = NULL, y = "Mean Fibrosis", title = NULL) +
#   scale_fill_manual(values = c("SHAM" = "grey", "ORAB" = "coral")) +
#   scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
#   theme_minimal()
# 
# 
# #ggsave(file.path(output_dir_plot, paste("fibrosis_barplot.png", sep = "")), plot = fibrosis_barplot, width = 8, height = 6)
# ggsave(file.path(output_dir_plot, paste("fibrosis_barplot.pdf", sep = "")), plot = fibrosis_barplot, width = 8, height = 6)
