# Create plot of percentage stressed cells
# Silje Wilhelmsen


# Visualization of stressed cardiomyocytes
##################################################################################
# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(writexl)
library(patchwork)
library(readr)

mouse_vcm_all_time_points_with_stress_status <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_vcm_all_time_points_with_stress_status.Rds")
output_dir_plot <- "C:/Users/siljeew/Master_project/snRNAseq/Plots"


# Plot of CM of total cell types 
##################################################################################
all_cells_count <- read_csv("C:/Users/siljeew/Master_project/snRNAseq/Data/pct_per_sample.csv")
# all_cells_count <- read_csv("C:/Users/siljeew/Master_project/snRNAseq/Data/pct_per_sample_25_05_12.csv")

View(all_cells_count)

sham_samples <- c('Sample1-TP75B', 'Sample3-TP52A', 'TP3-1', 'TP3-4', 'TP4-1', 
                  'TP5-1', 'TP9-1', 'TP9-5', 'TP11-5', 'TP17-2','TP13-1','TP17-3','TP21-3','TP22-1','TP22-2','TP23-2','TP23-4','TP26-5')
orab_samples <- c('TP14-2', 'TP11-1A', 'Sample2-TP511B', 'Sample4-TP62A', 'TP10-2', 
                'TP9-2', 'TP52-2', 'TP58-4', 'TP2-4', 'TP51-3','TP14-1','TP15-3','TP19-2', 'TP19-4','TP22-3','TP23-1','TP24-5','TP26-4')

# Define the conditions for each timepoint
timepoint_06_hours <- c('TP19-2', 'TP19-4','TP21-3','TP22-1','TP22-2','TP22-3')
timepoint_12_hours <- c('TP23-1', 'TP23-2','TP23-4','TP24-5','TP26-4','TP26-5')
timepoint_24_hours <- c('TP17-2', 'TP14-2','TP13-1','TP14-1','TP17-3','TP15-3')
timepoint_1_week <- c('TP2-4', 'TP3-1', 'TP3-4', 'TP4-1', 'TP52-2', 'TP58-4')
timepoint_3_days <- c('TP9-1', 'TP9-2', 'TP9-5', 'TP10-2', 'TP11-5', 'TP11-1A')
timepoint_3_weeks <- c('TP5-1','Sample1-TP75B', 'Sample3-TP52A','Sample2-TP511B', 'Sample4-TP62A','TP51-3')

# Create the new 'Timepoint' column and assign values based on the conditions
all_cells_count <- all_cells_count %>%
  mutate(Timepoint = NA, 
         Condition = NA)

all_cells_count$Timepoint <- ifelse(all_cells_count$orig.ident %in% timepoint_24_hours, '24 Hours',
                                            ifelse(all_cells_count$orig.ident %in% timepoint_1_week, '1 Week',
                                                   ifelse(all_cells_count$orig.ident %in% timepoint_3_days, '3 Days',
                                                          ifelse(all_cells_count$orig.ident %in% timepoint_06_hours, '6 Hours',
                                                                 ifelse(all_cells_count$orig.ident %in% timepoint_12_hours, '12 Hours',
                                                                        '3 Weeks')))))

head(all_cells_count)

# Create the new 'Condition' column and assign values based on the condition
all_cells_count$Condition <- ifelse(all_cells_count$orig.ident %in% sham_samples, 'SHAM', 'ORAB')

all_cells_count <- all_cells_count %>%
  group_by(Timepoint, Condition) %>%
  mutate(Percentage_cms_of_total = (vent_nuclei/total_nuclei) * 100)

head(all_cells_count)

all_cells_count <- all_cells_count %>%
  filter(Timepoint != "NA")

all_cells_count$Timepoint <- factor(all_cells_count$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
all_cells_count$Condition <- factor(all_cells_count$Condition, levels = c("SHAM", "ORAB"))


all_cells_count_plot <- ggplot(all_cells_count, aes(x = Timepoint, y = pct_vent, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(title = "Distribution of CMs of all celltypes",
       x = NULL,
       y = "CM of total nuclei (%)",
       fill = "Condition") +
  scale_y_continuous(breaks = seq(0, 100, by = 10), labels = paste0(seq(0, 100, by = 10), "%")) +
  scale_fill_manual(values = c("SHAM" = "grey22", "ORAB" = "coral")) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 26, color = "black", angle = 30, vjust = 0.8, margin = margin(t = 15)),
        axis.text.y = element_text(size = 26, colour = "black"),
        axis.title.y = element_text(size = 26, colour = "black", margin = margin(r = 15)),
        plot.title = element_text(size = 30, color = "black", margin = margin(b = 20, t = 10)),
        legend.title = element_text(size = 30, colour = "black"),
        legend.text = element_text(size = 30, color = "black"),
        axis.line.y.right = element_blank())

print(all_cells_count_plot)



# Plot of SHAM and ORAB of total count
##################################################################################
# Calculate nuclei count per group (ORAB and SHAM)
nuclei_count_data <- mouse_vcm_all_time_points_with_stress_status@meta.data %>%
  group_by(Timepoint, Group) %>%
  summarize(Cell_Count = n())

# Add total count per time point column
nuclei_count_data <- nuclei_count_data %>%
  group_by(Timepoint) %>%
  mutate(total_count = sum(Cell_Count))

# Calculate percentage nuclei of each group of total nuclei count
nuclei_count_data <- nuclei_count_data %>%
  mutate(percentage_of_total = (Cell_Count / total_count) * 100)

print(nuclei_count_data)  

nuclei_count_data$Timepoint <- factor(nuclei_count_data$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
nuclei_count_data$Group <- factor(nuclei_count_data$Group, levels = c("SHAM", "ORAB"))


# Create plot
nuclei_plot_condition <- ggplot(nuclei_count_data, aes(x = Timepoint, y = percentage_of_total, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(title = "Distribution of Nuclei",
       x = NULL,
       y = "Nuclei Count of Total Count (%)",
       fill = "Condition") +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
  scale_fill_manual(values = c("SHAM" = "grey22", "ORAB" = "coral")) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 26, color = "black", angle = 30, vjust = 0.8, margin = margin(t = 15)),
        axis.text.y = element_text(size = 26, colour = "black"),
        axis.title.y = element_text(size = 26, colour = "black", margin = margin(r = 15)),
        plot.title = element_text(size = 30, color = "black", margin = margin(b = 20, t = 10)),
        legend.title = element_text(size = 30, colour = "black"),
        legend.text = element_text(size = 30, color = "black"),
        axis.line.y.right = element_blank())


print(nuclei_plot_condition)

ggsave(file.path(output_dir_plot, paste("nuclei_plot_condition.png", sep = "")), plot = nuclei_plot_condition, width = 12, height = 10, dpi = 400)
ggsave(file.path(output_dir_plot, paste("nuclei_plot_condition.tiff", sep = "")), plot = nuclei_plot_condition, width = 12, height = 10, dpi = 400)





# Plot of SHAM and plot of ORAB combined
##################################################################################
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
  filter(Stress_Status != "SHAM - CM") %>%
  group_by(Timepoint, Stress_Status) %>%
  summarise(Cell_Count = n()) %>%
  group_by(Timepoint) %>%
  mutate(Total_Cells = sum(Cell_Count),
         Percentage_of_total_cells = (Cell_Count / Total_Cells) * 100)

print(percentage_data)

# percentage_data <- percentage_data %>%
#   mutate(Stress_Status = ifelse(Stress_Status == "SHAM - CM", "SHAM CM", Stress_Status))
# View(percentage_data)

# Set order of stress status
percentage_data$Stress_Status <- factor(percentage_data$Stress_Status,
                                        levels = c("Not stressed CM", "Stressed CM"))

percentage_data$Timepoint <- factor(percentage_data$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

cm_percentage_plot <- ggplot(percentage_data, aes(x = Timepoint, y = Percentage_of_total_cells, fill = Stress_Status)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  labs(title = "Distribution of stress status in ORAB",
       x = NULL,
       y = "Stress status (%)",
       fill = "Stress status") +
  scale_y_continuous(breaks = seq(0, 100, by = 10), labels = paste0(seq(0, 100, by = 10), "%")) +
  scale_fill_manual(values = c("Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 26, color = "black", angle = 30, vjust = 0.8, margin = margin(t = 15)),
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




# Combine two plots
############################################################################################################
all_cells_count_plot <- all_cells_count_plot +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))  # top, right, bottom, left margin

cm_percentage_plot <- cm_percentage_plot +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

stress_in_sham_plot <- stress_in_sham_plot +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
        
total_cm_distribution_plot <- (all_cells_count_plot | cm_percentage_plot) /
  (stress_in_sham_plot | plot_spacer())
total_cm_distribution_plot

ggsave(file.path(output_dir_plot, paste("cm_count_and_stress_distribution_plot.tiff", sep = "")), plot = total_cm_distribution_plot, width = 24, height = 20)

ggsave(file.path(output_dir_plot, paste("cm_count_and_stress_distribution_plot.png", sep = "")), plot = total_cm_distribution_plot, width = 24, height = 20)
ggsave(file.path(output_dir_plot, paste("cm_count_and_stress_distribution_plot.pdf", sep = "")), plot = total_cm_distribution_plot, width = 24, height = 20)


# write_xlsx(percentage_data_test, path = "C:/Users/siljeew/Master_project/snRNAseq/Data/cm_count.xlsx")
