# Boxplot 
# Input: dataset of small size (e.g. n=10 per group)
# Output: 
# Silje Wilhelmsen

# Preparation
###############################################################################
source("C:/Users/siljeew/Master_project/Phenotypic_data/Analysis/prepare_fibrosis_data.R") # Output: fibrosis_data


# Load libraries
library(openxlsx)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyverse)
library(writexl)
library(tidyr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(forcats)
library(cowplot)


# File path to the Excel file
fibrosis_data <- read.xlsx("C:/Users/siljeew/Master_project/Phenotypic_data/Data/Fibrosis.xlsx")
output_dir_plot <- "C:/Users/siljeew/Master_project/Phenotypic_data/Plots/"

unique(fibrosis_data$Timepoint)
fibrosis_data$timepoint <- trimws(fibrosis_data$Timepoint)


# Convert to factors and specify the desired order for timepoint
fibrosis_data$Timepoint <- factor(fibrosis_data$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
fibrosis_data$Condition <- factor(fibrosis_data$Condition, levels = c("SHAM", "ORAB"))
fibrosis_data$fibrosis <- as.numeric(fibrosis_data$fibrosis)

# Check levels of 'Timepoint'
levels(fibrosis_data$Timepoint)

 
# Boxplot
##########################################################################################
# Add significance data manually
significance_label <- data.frame(
  Timepoint = "3 Weeks",
  Condition = "SHAM",
  significance_label = "*",
  y_position = "5")
# View(significance_label)

significance_label$y_position <- as.numeric(significance_label$y_position)


fibrosis_boxplot <- ggplot(fibrosis_data, aes(x = Timepoint, y = fibrosis, fill = Condition)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8, lwd = 1.5, show.legend = FALSE) +
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), color = "black", size = 3, alpha = 0.8, show.legend = FALSE) +
  labs(title = NULL, 
       x = NULL, 
       y = "Fibrosis area (%)") +
  scale_fill_manual(values = c("SHAM" = "grey22", "ORAB" = "coral")) +
  scale_y_continuous(limits = c(0.2, 5.5)) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 32, color = "black", margin = margin(t = 15)), 
        axis.text.y = element_text(size = 32, colour = "black"), 
        axis.title.y = element_text(size = 50, margin = margin(r = 15)), # Change margins: r = right, b = bottom, t = top, l = left
        plot.title = element_text(size = 36, margin = margin(b = 20))) +
  geom_text(data = significance_label, aes(x = Timepoint, y = y_position, label = significance_label), 
            size = 32, vjust = 0.5, color = "black")

print(fibrosis_boxplot)


# Create a legend using a separate data frame
legend_df <- data.frame(
  Condition = c("SHAM", "ORAB"),
  color = c("grey22", "coral"),
  fibrosis = c(1, 1)
)

legend_df$Condition <- factor(legend_df$Condition, levels = c("SHAM", "ORAB"))

# Create a separate plot for legend
legend_plot <- ggplot(legend_df, aes(x = factor(1), y = fibrosis, fill = Condition)) +
  geom_bar(stat = "identity", width = 1, position = position_dodge()) +
  labs(fill = "Condition") +
  scale_fill_manual(values = c("SHAM" = "grey22", "ORAB" = "coral")) +
  theme_void() +
  theme(legend.title = element_text(size = 40),  
        legend.text = element_text(size = 40))

# Extract the legend using cowplot functions
legend <- get_legend(legend_plot)

# Combine the main plot with the manual legend
combined_plot <- plot_grid(fibrosis_boxplot, legend, ncol = 2, rel_widths = c(5, 1)) # The first plot is 4 times wider than the second

# Display the final plot
combined_plot

# ggsave(file.path(output_dir_plot, paste("fibrosis_boxplot.png", sep = "")), plot = combined_plot, width = 16, height = 12)

ggsave(file.path(output_dir_plot, paste("fibrosis_boxplot_small.png", sep = "")), plot = combined_plot, width = 15, height = 10)

# ggsave(file.path(output_dir_plot, paste("fibrosis_boxplot.pdf", sep = "")), plot = combined_plot, width = 8, height = 6, dpi = 400)




# Create graph
#######################################################
head(fibrosis_data)

fibrosis_data$`fibrosis` <- as.numeric(fibrosis_data$`fibrosis`)

# Create histogram
fibrosis_histogram <- fibrosis_data %>%
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

fibrosis_histogram



# Investigate data 
#####################################################################################
summary(fibrosis_data$fibrosis)
t.test(fibrosis_data)


# Barplot
##########################################################################################
# Prepare data for plot
data_for_plot <- kruskal_dunn_test_results %>%
  filter(group1_timepoint == group2_timepoint) %>%
  select(group1, group2, group1_timepoint, median_group1, median_group2 ,iqr_group1, iqr_group2, p_value_adjusted)

data_for_plot <- data_for_plot %>%
  rename(Timepoint = 'group1_timepoint',
         orab_median = 'median_group1',
         sham_median = 'median_group2',
         iqr_orab = 'iqr_group1',
         iqr_sham = 'iqr_group2') %>%
  pivot_longer(cols = c(orab_median, sham_median), 
               names_to = "Condition", 
               values_to = "median") %>%
  mutate(Condition = ifelse(Condition == "orab_median", "ORAB", "SHAM")) %>%
  mutate(iqr = ifelse(Condition == "ORAB", iqr_orab, iqr_sham))

data_for_plot <- data_for_plot %>%
  select(-group1, -group2, -iqr_orab, -iqr_sham)


data_for_plot <- data_for_plot %>%
  mutate(significance_label = ifelse(p_value_adjusted <0.05, "*", 
                                     ifelse(p_value_adjusted <0.01, "**",
                                             ifelse(p_value_adjusted <0.001, "***", 
                                                    "ns")))) 
print(data_for_plot)

significance_data <- data_for_plot %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = median + iqr + 0.1) 


# Set order of conditions
data_for_plot$Condition <- factor(data_for_plot$Condition, levels = c("SHAM", "ORAB"))
data_for_plot$Timepoint <- factor(data_for_plot$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))



# Create barplot with IQR-bar
#####################################################################
fibrosis_barplot <- ggplot(data_for_plot, aes(x = Timepoint, y = median, fill = Condition)) +
  geom_bar(stat = "identity", width = 0.9, position = position_dodge()) +
  geom_errorbar(aes(ymin = median - iqr, ymax = median + iqr), 
                width = 0.5, linewidth = 1,
                color = "black",
                position = position_dodge(0.9)) +
  labs(x = NULL, 
       y = "Fibrosis area (%)",
       title = "Fibrosis",
       fill = "Condition") +
  theme_classic() +
  scale_fill_manual(values = c("SHAM" = "grey22", "ORAB" = "coral")) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 18, color = "black", margin = margin(t = 15)), 
        axis.text.y = element_text(size = 16, colour = "black"), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 26, margin = margin(r = 15)), # Change margins: r = right, b = bottom, t = top, l = left
        plot.title = element_text(size = 26, margin = margin(b = 20)),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.line.y.right = element_blank()) +
  scale_y_continuous(limits = c(-0.1, 4))   # Add extra space above the bars
  
  #scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
  

# Display the plot
fibrosis_barplot

print(data_for_plot)

ggsave(file.path(output_dir_plot, paste("fibrosis_barplot.png", sep = "")), plot = fibrosis_barplot, width = 8, height = 6)
ggsave(file.path(output_dir_plot, paste("fibrosis_barplot.pdf", sep = "")), plot = fibrosis_barplot, width = 8, height = 6)
