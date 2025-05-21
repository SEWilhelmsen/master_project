# Cardiac remodelling plot for summary plot
# Silje Wilhelmsen

# Summary plot: 
file.edit("C:/Users/siljeew/Master_project/create_summary_plot.R")


library(Seurat)
library(writexl)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(readxl)
library(cowplot)
library(stringr)



# Fibrosis
###############################################################################
fibrosis_data <- read_excel("C:/Users/siljeew/Master_project/Phenotypic_data/Data/Fibrosis.xlsx")
head(fibrosis_data)

# Create data frame for plot 
data_remodelling <- data.frame(
  Timepoint = NA,
  Variable = NA,
  Percentage_to_sham = NA
)


# Calculate means for each timepoint and condition
mean_fibrosis <- fibrosis_data %>%
  group_by(timepoint, condition) %>%
  summarize(mean_fibrosis = mean(fibrosis, na.rm = TRUE), .groups = 'drop')

# Pivot the dataset to wide format
mean_fibrosis_wide <- mean_fibrosis %>%
  pivot_wider(names_from = condition, values_from = mean_fibrosis)

# Calculate the percentage of 'ORAB' mean to 'SHAM' mean
fibrosis_percentage <- mean_fibrosis_wide %>%
  mutate(percentage_to_sham = (ORAB / SHAM) * 100) %>%
  rename(mean_orab = ORAB,
         mean_sham = SHAM,
         Timepoint = timepoint)


# Populate the data frame with your percentage data 
data_remodelling <- data.frame(
  Timepoint = fibrosis_percentage$Timepoint,
  Variable = "Fibrosis",  # Set the variable for all rows
  Percentage_to_sham = fibrosis_percentage$percentage_to_sham
)

# Inspect the resulting data frame
print(data_remodelling)
data_remodelling$Timepoint <- trimws(as.character(data_remodelling$Timepoint))
data_remodelling$Timepoint <- factor(data_remodelling$Timepoint)
data_remodelling$Percentage_to_sham <- as.numeric(data_remodelling$Percentage_to_sham)




# LVW/BW
###################################################################################################
t_test_results <- read_csv("C:/Users/siljeew/Master_project/Phenotypic_data/Data/t_test_results.csv")
lvw_bw_data <- t_test_results %>%
  filter(variable == 'lvw_bw') 

View(lvw_bw_data)

# Assuming lvw_bw_data is already in the correct format, just select and rename necessary columns
lvw_bw_data <- lvw_bw_data %>%
  select(Timepoint, mean_sham, mean_orab) %>%  # Select required columns
  mutate(Variable = "LVW/BW",
         Percentage_to_sham = (mean_orab / mean_sham) * 100)

lvw_bw_data <- lvw_bw_data %>% 
  select(-mean_sham, -mean_orab)

# Merge into data frame
data_remodelling <- bind_rows(data_remodelling, lvw_bw_data)
print(data_remodelling)

# Set as factor for correct order
data_remodelling$Timepoint <- factor(data_remodelling$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
data_remodelling$Variable <- factor(data_remodelling$Variable)
data_remodelling$Percentage_to_sham <- as.numeric(data_remodelling$Percentage_to_sham)
print(data_remodelling)






# Nuclei count
####################################################################################
mouse_vcm_all_time_points_with_stress_status <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_vcm_all_time_points_with_stress_status.Rds")
output_dir_plot <- "C:/Users/siljeew/Master_project/snRNAseq/Plots"


nuclei_count <- mouse_vcm_all_time_points_with_stress_status@meta.data %>%
  group_by(Timepoint, Stress_Status) %>%
  summarize(Cell_Count = n()) %>%
  group_by(Timepoint, Stress_Status)

# Remove line in SHAM
nuclei_count <- nuclei_count %>%
  mutate(Stress_Status = ifelse(Stress_Status == "SHAM - CM", "SHAM CM", Stress_Status))


nuclei_count <- nuclei_count %>%
  group_by(Timepoint) %>%
  mutate(cell_count_sham = Cell_Count[Stress_Status == "SHAM CM"])

nuclei_count <- nuclei_count %>%
  mutate(percentage_to_sham = (Cell_Count / cell_count_sham) * 100)
# print(nuclei_count)

# Percentage stressed of all ORAB CMs
orab_nuclei <- nuclei_count %>%
  filter(Stress_Status != "SHAM CM") %>%
  dplyr::select(Timepoint, Stress_Status, Cell_Count)

orab_nuclei <- orab_nuclei %>%
  group_by(Timepoint) %>%
  mutate(total_orab = sum(Cell_Count))
  

orab_nuclei <- orab_nuclei %>%
  mutate(percentage_of_total = (Cell_Count/total_orab) * 100) 

nuclei_data <- orab_nuclei %>%
  filter(Stress_Status == "Stressed CM") %>%
  select(Timepoint, Stress_Status, percentage_of_total) %>%
  rename(Variable = Stress_Status,
         Percentage_to_sham = percentage_of_total)
  

nuclei_data <- nuclei_data %>%
  mutate(Variable = ifelse(Variable == "Stressed CM", "% of CMs that are stressed", Variable))
print(nuclei_data)

data_remodelling <- data_remodelling %>%
  filter(Variable != "NA")

# data_remodelling <- data_remodelling %>%
#   filter(Variable != "Stressed CM",
#          Variable != "Not stressed CM",
#          Variable != "% CM of that are stressed")

# Merge into data frame
data_remodelling <- bind_rows(data_remodelling, nuclei_data)
print(data_remodelling)

# Set as factor for correct order
data_remodelling$Timepoint <- factor(data_remodelling$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
data_remodelling$Variable <- factor(data_remodelling$Variable, levels = c("Fibrosis", "LVW/BW", "% of CMs that are stressed"))
data_remodelling$Percentage_to_sham <- as.numeric(data_remodelling$Percentage_to_sham)
print(data_remodelling)

# plot <- ggplot(orab_nuclei, aes(x = Timepoint, y = percentage_of_total, fill = Stress_Status)) +
#   geom_bar(stat = "identity", position = "stack") +
#   theme_classic() +
#   labs(title = "Nuclei Percentage for Summary Figure",
#        x = NULL,
#        y = "Nuclei Count compared to SHAM (%)",
#        fill = "Stress Status") +
#   scale_fill_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
#   theme(axis.text.x = element_text(hjust = 0.5, size = 26, color = "black", angle = 30, vjust = 0.8, margin = margin(t = 15)),
#         axis.text.y = element_text(size = 26, colour = "black"),
#         axis.title.y = element_text(size = 26, colour = "black", margin = margin(r = 15)),
#         plot.title = element_text(size = 30, color = "black", margin = margin(b = 20, t = 10)),
#         legend.title = element_text(size = 30, colour = "black"),
#         legend.text = element_text(size = 30, color = "black"),
#         axis.line.y.right = element_blank())
# 
# 
# print(plot)


# # Create plot
# label_data <- data_remodelling %>%
#   filter(Timepoint == "3 Weeks") %>%
#   mutate(Percentage_to_sham = Percentage_to_sham + 1)  # Adjust the y position above the point

# Create the plot with line and point representation
remodelling_plot <- ggplot(data_remodelling, aes(x = Timepoint, y = Percentage_to_sham, color = Variable, group = Variable)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  geom_hline(yintercept = 100, linetype = "dashed", color = "black") +  # Add horizontal line at y = 100
  # geom_text(data = label_data, aes(x = Timepoint, y = Percentage_to_sham, label = Variable), 
  #           color = "black", size = 6, vjust = -0.5, hjust = 0.5) +  # Position the text at the appropriate y position
  labs(title = "Cardiac remodelling and morphology", 
       x = "Time", 
       y = "Values compared to SHAM (%)") +
  scale_color_manual(
    values = c("Fibrosis" = "darkcyan", 
               "LVW/BW" = "brown", 
               "% of CMs that are stressed" = "olivedrab4"), 
    labels = c("Fibrosis", "LVW/BW", "% of CMs \nthat are stressed")  # \n gives a lineshift
  ) +
  scale_y_continuous(breaks = seq(0, 240, by = 30), limits = c(0, 240))  +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 26, color = "black", margin = margin(t = 15)),
        axis.text.y = element_text(size = 26, colour = "black"), 
        axis.title.y = element_text(size = 26, margin = margin(r = 15)),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 26, margin = margin(b = 20)),
        legend.title = element_text(size = 26, colour = "black"),
        legend.text = element_text(size = 26, colour = "black"),
        #legend.position = "none",
        plot.margin = margin(t = 20, r = 30, b = 10, l = 10))

# Display the plot
print(remodelling_plot)


