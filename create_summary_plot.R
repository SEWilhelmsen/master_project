# Master summary figure 
# Silje Wilhelmsen

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(readxl)



# Cardiac remodelling plot (fibrosis, cell size?)
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

# Load cell size data
cell_size_data_filtered <- read_excel("C:/Users/siljeew/Master_project/Phenotypic_data/Data/cell_size_data_filtered.xlsx")

# Calculate means for each timepoint and condition
mean_cell_size <- cell_size_data_filtered %>%
  group_by(Timepoint, Condition) %>%
  summarize(mean_area = mean(area_um2, na.rm = TRUE), .groups = 'drop')

# Pivot the dataset to wide format
mean_cell_size_wide <- mean_cell_size %>%
  pivot_wider(names_from = Condition, values_from = mean_area)

# Calculate the percentage of 'ORAB' mean to 'SHAM' mean
area_percentage <- mean_cell_size_wide %>%
  mutate(percentage_to_sham = (ORAB / SHAM) * 100) %>%
  rename(mean_orab = ORAB,
         mean_sham = SHAM,
         Timepoint = Timepoint)

# Create data frame with your percentage data 
data_area <- data.frame(
  Timepoint = area_percentage$Timepoint,
  Variable = "Cell Area",  # Set the variable for all rows
  Percentage_to_sham = area_percentage$percentage_to_sham
)

# Combine both data frames by adding new rows
data_remodelling <- rbind(data_remodelling, data_area)
print(data_remodelling)


# LVW/BW
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



# Create plot
label_data <- data_remodelling %>%
  filter(Timepoint == "3 Weeks") %>%
  mutate(Percentage_to_sham = Percentage_to_sham + 1)  # Adjust the y position above the point

# Create the plot with line and point representation
remodelling_plot <- ggplot(data_remodelling, aes(x = Timepoint, y = Percentage_to_sham, color = Variable, group = Variable)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  geom_hline(yintercept = 100, linetype = "dashed", color = "black") +  # Add horizontal line at y = 100
  # geom_text(data = label_data, aes(x = Timepoint, y = Percentage_to_sham, label = Variable), 
  #           color = "black", size = 6, vjust = -0.5, hjust = 0.5) +  # Position the text at the appropriate y position
  labs(title = "Cardiac Remodelling and Morphology", 
       x = "Time", 
       y = "Values compared to SHAM (%)") +
  scale_color_manual(values = c("Fibrosis" = "darkcyan", "Cell Area" = "goldenrod2", "LVW/BW" = "brown")) +
  scale_y_continuous(breaks = seq(80, 240, by = 20), limits = c(80, 240))  +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 26, color = "black", margin = margin(t = 15)),
        axis.text.y = element_text(size = 26, colour = "black"), 
        axis.title.y = element_text(size = 26, margin = margin(r = 15)),
        axis.title.x = element_text(hjust = 1, vjust = 1, size = 26, margin = margin(t = 15)),
        plot.title = element_text(size = 26, margin = margin(b = 20)),
        legend.title = element_text(size = 26, colour = "black"),
        legend.text = element_text(size = 26, colour = "black"),
        #legend.position = "none",
        plot.margin = margin(t = 10, r = 100, b = 10, l = 10))

# Display the plot
print(remodelling_plot)




# Transcription of energy metabolic genes in CM plot
##############################################################################
# Function to populate data_transcription with the results of a specific process.
populate_data_transcription <- function(process_name) {

  t_test_go_process_results <- read_excel("C:/Users/siljeew/Master_project/snRNAseq/Data/t_test_go_process_results.xlsx")
  
  # Filter for the specified process and remove NA and unwanted time point
  process_data <- t_test_go_process_results %>%
    filter(process == process_name & mean_orab != "NA") %>%
    filter(time_point != "1 week") # Filter out unwanted value
  
  # Select and calculate Percentage to sham
  process_data <- process_data %>%
    select(time_point, mean_sham, mean_orab) %>%
    mutate(Percentage_to_sham = (mean_orab / mean_sham) * 100) %>%
    select(-mean_sham, -mean_orab) # Drop unnecessary columns
  
  # Create a data frame with populated data
  process_df <- data.frame(
    Timepoint = process_data$time_point,
    Process = process_name,
    Percentage_to_sham = process_data$Percentage_to_sham
  )
  
  return(process_df)
}

# Initialize data_transcription data frame
data_transcription <- data.frame(
  Timepoint = NA,
  Process = NA,
  Percentage_to_sham = NA
)

# TCA
tca_data <- populate_data_transcription("tca")
data_transcription <- bind_rows(data_transcription, tca_data)

# B-Oxidation
b_oxidation_data <- populate_data_transcription("b_oxidation")
data_transcription <- bind_rows(data_transcription, b_oxidation_data)

# Glycolysis
glycolysis_data <- populate_data_transcription("glycolysis")
data_transcription <- bind_rows(data_transcription, glycolysis_data)

# Ketones
ketones_data <- populate_data_transcription("ketone_catabolism")
data_transcription <- bind_rows(data_transcription, ketones_data)

# Print the final data_transcription
print(data_transcription)

data_transcription <- data_transcription %>% 
  mutate(Process = case_when(
    Process == "tca" ~"TCA Cycle", 
    Process == "glycolysis" ~"Glycolysis",
    Process == "b_oxidation" ~"B-Oxidation",
    Process == "ketone_catabolism" ~"Ketone Catabolism",
    TRUE ~Process # The rest of the process names are kept unchanged
  )) %>%
  filter(!is.na(Percentage_to_sham))

data_transcription$Timepoint <- factor(data_transcription$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
data_transcription$Process <- factor(data_transcription$Process)
data_transcription$Percentage_to_sham <- as.numeric(data_transcription$Percentage_to_sham)


# Create plot
# label_data <- data_transcription %>%
#   filter(Timepoint == "3 Weeks") %>%
#   mutate(Percentage_to_sham = Percentage_to_sham + 1)  # Adjust the y position above the point


transcription_plot <- ggplot(data_transcription, aes(x = Timepoint, y = Percentage_to_sham, color = Process, group = Process)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  geom_hline(yintercept = 100, linetype = "longdash", linewidth = 1, color = "grey22") +  # Add horizontal line at y = 100
  labs(title = "Transcription of Energy Metabolic Genes in Cardiomyocytes", 
       x = "Time", 
       y = "Values compared to SHAM (%)") +
  scale_color_manual(values = c("TCA Cycle" = "darkcyan", "Glycolysis" = "goldenrod2", 
                                "B-Oxidation" = "brown", "Ketone Catabolism" = "chocolate1")) +
  scale_y_continuous(breaks = seq(70, 130, by = 10), limits = c(70, 130))  +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 26, color = "black", margin = margin(t = 15)),
        axis.text.y = element_text(size = 26, colour = "black"), 
        axis.title.y = element_text(size = 26, margin = margin(r = 15)),
        axis.title.x = element_text(hjust = 1, vjust = 1, size = 26, margin = margin(t = 15)),
        plot.title = element_text(size = 26, margin = margin(b = 20)),
        legend.title = element_text(size = 26, colour = "black"),
        legend.text = element_text(size = 26, colour = "black"),
        #legend.position = "none",
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

transcription_plot 



# Activation of signalling pathways 
###############################################################################
# Load data
signal_diff_data <- read_excel("C:/Users/siljeew/Master_project/snRNAseq/Data/tukey_results_progeny_by_condition.xlsx")
View(signal_diff_data)

# Define activity of SHAM as 1 for relative comparison
activity_sham_value <- 1

signal_diff_data <- signal_diff_data %>%
  mutate(
    Activity_orab = activity_sham_value + diff,
    Percentage_to_sham = (Activity_orab / activity_sham_value) * 100
  )
print(signal_diff_data)

signal_diff_data$Timepoint <- factor(signal_diff_data$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Filter out TNFa because it is aligned with NFkB
signal_diff_data <- signal_diff_data %>%
  filter(Pathway != "TNFa",
         Pathway != "Androgen",)

signal_plot <- ggplot(signal_diff_data, aes(x = Timepoint, y = Percentage_to_sham, color = Pathway, group = Pathway)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  geom_hline(yintercept = 100, linetype = "longdash", linewidth = 1, color = "grey22") +  # Add horizontal line at y = 100
  labs(title = "Activity of Signaling Pathways in Cardiomyocytes", 
       x = "Time", 
       y = "Values compared to SHAM (%)") +
  scale_color_manual(values = c("p53" = "darkcyan", "EGFR" = "goldenrod2", "NFkB" = "brown", "JAK-STAT" = "olivedrab4")) +
  scale_y_continuous(breaks = seq(80, 115, by = 5), limits = c(80, 115))  +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 26, color = "black", margin = margin(t = 15)),
        axis.text.y = element_text(size = 26, colour = "black"), 
        axis.title.y = element_text(size = 26, margin = margin(r = 15)),
        axis.title.x = element_text(hjust = 1, vjust = 1, size = 26, margin = margin(t = 15)),
        plot.title = element_text(size = 26, margin = margin(b = 20)),
        legend.title = element_text(size = 26, colour = "black"),
        legend.text = element_text(size = 26, colour = "black"),
        plot.margin = margin(t = 10, r = 100, b = 10, l = 10))

signal_plot


# Combine plot 
###############################################################################
combined_plot <- plot_grid(remodelling_plot, transcription_plot, signal_plot, ncol = 1, nrow = 3)
combined_plot

ggsave("C:/Users/siljeew/Master_project/summary_plot.png", plot = combined_plot, height = 20, width = 15)
ggsave("C:/Users/siljeew/Master_project/summary_plot.tiff", plot = combined_plot, height = 20, width = 15)