# Activation of signalling pathways plot for summary plot
# Silje Wilhelmsen

# Summary plot: 
file.edit("C:/Users/siljeew/Master_project/create_summary_plot.R")

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(readxl)
library(cowplot)
library(stringr)


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

signal_diff_data <- signal_diff_data %>%
  mutate(Pathway = str_replace(Pathway, "NFkB", "NF\u03BAB, TNF\u03B1")) %>%
  mutate(Pathway = str_replace(Pathway, "p53", "p53, Androgen"))

# Filter out TNFa because it is aligned with NFkB
signal_diff_data <- signal_diff_data %>%
  filter(Pathway != "TNFa",
         Pathway != "Androgen")

print(signal_diff_data)

signal_diff_data$Timepoint <- factor(signal_diff_data$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))



signal_plot <- ggplot(signal_diff_data, aes(x = Timepoint, y = Percentage_to_sham, color = Pathway, group = Pathway)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  geom_hline(yintercept = 100, linetype = "longdash", linewidth = 1, color = "grey22") +  # Add horizontal line at y = 100
  labs(title = "Activity of Signaling Pathways in Cardiomyocytes", 
       x = "Time", 
       y = "Values compared to SHAM (%)") +
  scale_color_manual(values = c("p53, Androgen" = "darkcyan", 
                                "EGFR" = "goldenrod2", 
                                "NF\u03BAB, TNF\u03B1" = "brown", 
                                "JAK-STAT" = "olivedrab4")) +
  scale_y_continuous(breaks = seq(80, 115, by = 5), limits = c(80, 115))  +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 26, color = "black", margin = margin(t = 15)),
        axis.text.y = element_text(size = 26, colour = "black"), 
        axis.title.y = element_text(size = 26, margin = margin(r = 15)),
        axis.title.x = element_text(hjust = 1, vjust = 1, size = 26, margin = margin(t = 15)),
        plot.title = element_text(size = 26, margin = margin(b = 20)),
        legend.title = element_text(size = 26, colour = "black"),
        legend.text = element_text(size = 26, colour = "black"),
        plot.margin = margin(t = 10, r = 50, b = 10, l = 10))

signal_plot