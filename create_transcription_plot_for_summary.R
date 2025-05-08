# Transcription of energy metabolic genes in summary plot
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
                                "B-Oxidation" = "brown", "Ketone Catabolism" = "olivedrab4")) +
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
