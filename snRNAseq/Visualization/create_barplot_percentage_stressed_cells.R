# Create plot of percentage stressed cells
# Silje Wilhelmsen


# Visualization of stressed cardiomyocytes
##################################################################################
# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(writexl)

mouse_vcm_all_time_points_with_stress_status <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_vcm_all_time_points_with_stress_status.Rds")
output_dir_plot <- "C:/Users/siljeew/Master_project/snRNAseq/Plots"


# Plot of SHAM and plot of ORAB combined
##################################################################################
# nuclei_sham <- mouse_vcm_all_time_points_with_stress_status@meta.data %>%
#   group_by(Timepoint, Stress_Status) %>%
#   filter(Stress_Status == "SHAM - CM") %>%
#   summarize(Cell_Count = n()) %>%
#   group_by(Timepoint)
# 
# # Rename for better layout in plot 
# nuclei_sham <- nuclei_sham %>%
#   mutate(Stress_Status = ifelse(Stress_Status == "SHAM - CM", "SHAM CM", Stress_Status))
# 
# nuclei_sham$Timepoint <- factor(nuclei_sham$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
# 
# # Create plot for SHAM
# 
# 
# nuclei_orab <- mouse_vcm_all_time_points_with_stress_status@meta.data %>%
#   group_by(Timepoint, Stress_Status) %>%
#   filter(Stress_Status != "SHAM - CM") %>%
#   summarize(Cell_Count = n()) %>%
#   group_by(Timepoint)
# nuclei_orab$Timepoint <- factor(nuclei_orab$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
# nuclei_orab$Stress_Status <- factor(nuclei_orab$Stress_Status)


# Prepare data
nuclei_count <- mouse_vcm_all_time_points_with_stress_status@meta.data %>%
  group_by(Timepoint, Stress_Status) %>%
  summarize(Cell_Count = n()) %>%
  group_by(Timepoint, Stress_Status)

# Remove line in SHAM
nuclei_count <- nuclei_count %>%
  mutate(Stress_Status = ifelse(Stress_Status == "SHAM - CM", "SHAM CM", Stress_Status))

nuclei_count$Timepoint <- factor(nuclei_count$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
nuclei_count$Stress_Status <- factor(nuclei_count$Stress_Status, levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))

# Create plot
nuclei_count_plot <- ggplot(nuclei_count, aes(x = Timepoint, y = Cell_Count, fill = Stress_Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(title = "Distribution of Stress Status",
       x = NULL,
       y = "Nuclei Count",
       fill = "Stress Status") +
  scale_y_continuous(breaks = seq(0, 6000, by = 1000), limits = c(0, 6000)) +
  scale_fill_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 26, color = "black", angle = 30, vjust = 0.8, margin = margin(t = 15)),
        axis.text.y = element_text(size = 26, colour = "black"),
        axis.title.y = element_text(size = 26, colour = "black", margin = margin(r = 15)),
        plot.title = element_text(size = 30, color = "black", margin = margin(b = 20, t = 10)),
        legend.title = element_text(size = 30, colour = "black"),
        legend.text = element_text(size = 30, color = "black"),
        axis.line.y.right = element_blank())


print(nuclei_count_plot)

ggsave(file.path(output_dir_plot, paste("nuclei_count_plot.png", sep = "")), plot = nuclei_count_plot, width = 12, height = 10, dpi = 400)
ggsave(file.path(output_dir_plot, paste("nuclei_count_plot.tiff", sep = "")), plot = nuclei_count_plot, width = 12, height = 10, dpi = 400)



# Percentage plot
##################################################################################
percentage_data <- mouse_vcm_all_time_points_with_stress_status@meta.data %>%
  group_by(Timepoint, Stress_Status) %>%
  summarise(Cell_Count = n()) %>%
  group_by(Timepoint) %>%
  mutate(Total_Cells = sum(Cell_Count),
         Percentage_of_total_cells = (Cell_Count / Total_Cells) * 100)

print(percentage_data)

percentage_data <- percentage_data %>%
  mutate(Stress_Status = ifelse(Stress_Status == "SHAM - CM", "SHAM CM", Stress_Status))
View(percentage_data)

# Set order of stress status
percentage_data$Stress_Status <- factor(percentage_data$Stress_Status,
                                        levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))

percentage_data$Timepoint <- factor(percentage_data$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

cm_percentage_plot <- ggplot(percentage_data, aes(x = Timepoint, y = Percentage_of_total_cells, fill = Stress_Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(title = "Distribution of Stress Status",
       x = NULL,
       y = "Stress Status (%)",
       fill = "Stress Status") +
  scale_y_continuous(breaks = seq(0, 100, by = 10), labels = paste0(seq(0, 100, by = 10), "%")) +
  scale_fill_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 20, color = "black", margin = margin(t = 15)),
        axis.text.y = element_text(size = 26, colour = "black"),
        axis.title.y = element_text(size = 26, colour = "black", margin = margin(r = 15)),
        plot.title = element_text(size = 30, color = "black", margin = margin(b = 20, t = 10)),
        legend.title = element_text(size = 30, colour = "black"),
        legend.text = element_text(size = 30, color = "black"),
        axis.line.y.right = element_blank()) 

print(cm_percentage_plot)

ggsave(file.path(output_dir_plot, paste("cm_percentage_plot.png", sep = "")), plot = cm_percentage_plot, width = 12, height = 10, dpi = 400)
ggsave(file.path(output_dir_plot, paste("cm_percentage_plot.pdf", sep = "")), plot = cm_percentage_plot, width = 12, height = 10, dpi = 400)


# Plot with counts labels
########################################################################################
total_counts <- percentage_data %>%
  filter(Stress_Status == "SHAM CM")  
print(total_counts)

total_counts$Timepoint <- factor(total_counts$Timepoint)
percentage_data$Timepoint <- factor(percentage_data$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
percentage_data$Stress_Status <- factor(percentage_data$Stress_Status,
                                        levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))
percentage_data$Percentage_of_total_cells <- as.numeric(percentage_data$Percentage_of_total_cells)



# Use position = "stack" to present on top of each other
cm_percentage_plot <- ggplot(percentage_data, aes(x = Timepoint, y = Percentage_of_total_cells, fill = Stress_Status)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  labs(title = "Distribution of Stress Status",
       x = NULL,
       y = "Stress Status (%)",
       fill = "Stress Status") +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  scale_fill_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 26, color = "black", angle = 30, vjust = 0.8, margin = margin(t = 15)),
        axis.text.y = element_text(size = 26, colour = "black"),
        axis.title.y = element_text(size = 26, colour = "black", margin = margin(r = 15)),
        plot.title = element_text(size = 30, color = "black", margin = margin(b = 20, t = 10)),
        legend.title = element_text(size = 30, colour = "black"),
        legend.text = element_text(size = 30, color = "black"),
        axis.line.y.right = element_blank()) +
  geom_text(aes(label = Cell_Count), size = 8, position = position_stack(.5)) +
  geom_text(data = total_counts, aes(x = Timepoint, label = Total_Cells, y = 102), size = 8, 
            position = position_dodge(0.9), vjust = 0)

print(cm_percentage_plot)

ggsave(file.path(output_dir_plot, paste("cm_percentage_plot_stack.png", sep = "")), plot = cm_percentage_plot, width = 12, height = 10, dpi = 400)
ggsave(file.path(output_dir_plot, paste("cm_percentage_plot_stack.pdf", sep = "")), plot = cm_percentage_plot, width = 12, height = 10, dpi = 400)



# Investigate development for each Stress_Status compared to itself at 6 Hours
######################################################################################
# Calculate percentage of Cell_Count compared to 6 Hours
percentage_data_test <- percentage_data %>% 
  group_by(Stress_Status) %>% 
  mutate(baseline_value = Cell_Count[Timepoint == "6 Hours"]) %>%  
  mutate(percentage_of_6h = (Cell_Count / baseline_value) * 100) %>% 
  ungroup() %>%  # Remove grouping
  select(-baseline_value) 

# Print the modified data frame
print(percentage_data_test)



# Create plot
change_from_6h_plot <- ggplot(percentage_data_test, aes(x = Timepoint, y = percentage_of_6h, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 4) +
  labs(title = "Nuclei Count Compared to 6 Hours", 
       x = NULL,
       y = "Nuclei Count (%)",
       color = "Stress Status") +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  geom_hline(yintercept = 100) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 26, color = "black", angle = 30, vjust = 0.8, margin = margin(t = 15)),
        axis.text.y = element_text(size = 26, colour = "black"),
        axis.title.y = element_text(size = 26, colour = "black", margin = margin(r = 15)),
        plot.title = element_text(size = 30, color = "black", margin = margin(b = 20, t = 10)),
        legend.title = element_text(size = 30, colour = "black"),
        legend.text = element_text(size = 30, color = "black"),
        legend.position = "none", # Remove legends for combined plot
        axis.line.y.right = element_blank())

change_from_6h_plot


ggsave(file.path(output_dir_plot, paste("cm_percentage_plot_from_6h.png", sep = "")), plot = change_from_6h_plot, width = 12, height = 10, dpi = 400)
ggsave(file.path(output_dir_plot, paste("cm_percentage_plot_from_6h.pdf", sep = "")), plot = change_from_6h_plot, width = 12, height = 10, dpi = 400)


change_from_6h_plot <- change_from_6h_plot +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))  # top, right, bottom, left margin

cm_percentage_plot <- cm_percentage_plot +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

        
total_cm_distribution_plot <- (change_from_6h_plot | cm_percentage_plot)
total_cm_distribution_plot

ggsave(file.path(output_dir_plot, paste("total_cm_distribution_plot.png", sep = "")), plot = total_cm_distribution_plot, width = 24, height = 10, dpi = 400)
ggsave(file.path(output_dir_plot, paste("total_cm_distribution_plot.pdf", sep = "")), plot = total_cm_distribution_plot, width = 24, height = 10, dpi = 400)


# write_xlsx(percentage_data_test, path = "C:/Users/siljeew/Master_project/snRNAseq/Data/cm_count.xlsx")
