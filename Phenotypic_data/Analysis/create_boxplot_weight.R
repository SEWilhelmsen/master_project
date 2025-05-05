# Create boxplot for weight data 
# Silje Wilhelmsen

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
# processed_data <- read.xlsx("C:/Users/siljeew/Master_project/Phenotypic_data/Data/Fibrosis.xlsx")
output_dir_plot <- "C:/Users/siljeew/Master_project/Phenotypic_data/Plots/"

data <- processed_data

# Convert to factors and specify the desired order for timepoint
processed_data$Timepoint <- factor(processed_data$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
processed_data$Condition <- factor(processed_data$Condition, levels = c("SHAM", "ORAB"))
processed_data$lvw_bw <- as.numeric(processed_data$lvw_bw)
processed_data$`lw` <- as.numeric(processed_data$`lw`)

# Verify correct levels of 'Timepoint'
levels(processed_data$Timepoint)


# Boxplot LVW/BW
##########################################################################################
# Add significance data manually
significance_label <- data.frame(
  Timepoint = c("1 Day", "3 Days", "1 Week", "3 Weeks"),
  Condition = c("SHAM"),
  significance_label = c("**", "**", "***", "***"),
  y_position = c("5"))

significance_label$Condition <- factor(significance_label$Condition)
significance_label$Timepoint <- factor(significance_label$Timepoint)
significance_label$y_position <- as.numeric(significance_label$y_position)

lvw_bw_boxplot <- ggplot(processed_data, aes(x = Timepoint, y = lvw_bw, fill = Condition)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8, lwd = 1.5, show.legend = FALSE) +
  geom_jitter(data = filter(processed_data, Used_already_for == "Multiome seq"),
              aes(x = Timepoint, y = lvw_bw, color = Used_already_for, group = Condition),
              size = 6, shape = 21, fill = "red",
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
              show.legend = FALSE) +
  labs(title = NULL, 
       x = NULL, 
       y = "LVW/BW") +
  scale_fill_manual(values = c("SHAM" = "grey22", "ORAB" = "coral")) +
  scale_color_manual(name = "", values = c("Multiome seq" = "red"), labels = c("Sample used for RNA sequencing")) +
  scale_y_continuous(limits = c(2.5, 5)) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 32, color = "black", angle = 30, margin = margin(t = 15)), 
        axis.text.y = element_text(size = 32, colour = "black"), 
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 32, margin = margin(r = 15)), 
        plot.title = element_text(size = 32, margin = margin(b = 20)),
        axis.line.y.right = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = NA), title = "Condition"),  
         color = guide_legend(override.aes = list(shape = 21, fill = "red"), title = "Sample")) + 
  geom_text(data = significance_label, aes(x = Timepoint, y = y_position, label = significance_label),
            size = 24, vjust = 0.5, color = "black")

print(lvw_bw_boxplot)


# Boxplot LW
##########################################################################################
lw_boxplot <- ggplot(processed_data, aes(x = Timepoint, y = lw, fill = Condition)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.8, lwd = 1.5, show.legend = FALSE) +
  geom_jitter(data = filter(processed_data, Used_already_for == "Multiome seq"),
              aes(x = Timepoint, y = lw, color = Used_already_for, group = Condition),
              size = 5, shape = 21, fill = "red",
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              show.legend = FALSE) +
  labs(title = NULL, 
       x = NULL, 
       y = "LW (mg)") +
  scale_fill_manual(values = c("SHAM" = "grey22", "ORAB" = "coral")) +
  scale_color_manual(name = "", values = c("Multiome seq" = "red"), labels = c("Sample used for RNA sequencing")) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 32, color = "black", angle = 30, margin = margin(t = 15)), 
        axis.text.y = element_text(size = 32, colour = "black"), 
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 32, margin = margin(r = 15)), 
        plot.title = element_blank(),
        axis.line.y.right = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = NA), title = "Condition"),  
         color = guide_legend(override.aes = list(shape = 21, fill = "red"), title = "Sample"))
  

print(lw_boxplot)



# Create legends
####################################################################################
# Create a data frame for the legend with condition bars
legend_df <- data.frame(
  Condition = c("SHAM", "ORAB"),
  color = c("grey22", "coral"),
  value = c(0, 0)  # Bar height for condition
)
# Set correct order of conditions
legend_df$Condition <- factor(legend_df$Condition, levels = c("SHAM", "ORAB")) 

legend_plot <- ggplot() +
  geom_bar(data = legend_df, aes(x = Condition, y = value, fill = Condition), stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("SHAM" = "grey22", "ORAB" = "coral")) +
  geom_point(aes(x = 0, y = 0, color = "Samples used for\n RNA sequencing"), 
             shape = 21, size = 12, fill = "red") + # use \n for line shift
  labs(fill = "Condition", 
       color = NULL) + # color = NULL removes legend title of red points
  theme_void() +
  theme(legend.title = element_text(size = 42, colour = "black"),
        legend.text = element_text(size = 42, colour = "black")) +
  guides(fill = guide_legend(order = 1),
         color = guide_legend(order = 2))

# legend_plot

# Extract the legends 
legend <- get_legend(legend_plot)  



# Combine plots
################################################################################
# Set size and margins for plot
lvw_bw_boxplot <- lvw_bw_boxplot + theme(plot.margin = margin(t = 0, r = 10, b = 0, l = 10))
lw_boxplot <- lw_boxplot + theme(plot.margin = margin(t = 0, r = 10, b = 0, l = 10))

combined_plot <- ggdraw() + 
  draw_plot(lvw_bw_boxplot, x = 0, y = 0, width = 0.4, height = 1) +  # First plot
  draw_plot(lw_boxplot, x = 0.4, y = 0, width = 0.4, height = 1) +  # Second plot
  draw_plot(legend, x = 0.8, y = 0, width = 0.2, height = 1)  # Legend
  
# Print or save the combined plot
print(combined_plot)


# Combine LVW/BW plot, LW plot and legend plot
# combined_plot <- plot_grid(lvw_bw_boxplot, lw_boxplot,legend, ncol = 3, rel_widths = c(2,2,1.5))
# combined_plot <- plot_grid(lvw_bw_boxplot, lw_boxplot,legend, ncol = 3)
print(combined_plot)


ggsave(file.path(output_dir_plot, "weight_boxplot.tiff"), plot = combined_plot, width = 28, height = 10, dpi = 400)
ggsave(file.path(output_dir_plot, "weight_boxplot.png"), plot = combined_plot, width = 28, height = 10, dpi = 400)

