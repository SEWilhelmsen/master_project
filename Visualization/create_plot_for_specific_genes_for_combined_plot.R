# Create plot of specific genes to combine with plot of go_process
# Silje Wilhelmsen

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# Combine the data processing and plotting into a single function for simplicity

create_go_genes_plot <- function(data_avg_log2fc, go_process_of_interest_genes) {  
  
  # # Import data 
  # data_avg_log2fc <- read_csv2(file_path_to_data_avg_log2fc)
  #   
  # # Ensure gene names are uppercase
  # data_avg_log2fc$gene <- toupper(data_avg_log2fc$gene)
  
  # Subset the data to include only genes in go_process_of_interest_genes
  data_avg_log2fc_genes_of_interest <- data_avg_log2fc %>%
    filter(gene %in% go_process_of_interest_genes)
  
  # Reshape from wide to long format
  data_avg_log2fc_genes_of_interest_go_process_long <- data_avg_log2fc_genes_of_interest %>%
    pivot_longer(cols = -gene, names_to = "time_point", values_to = "avg_log2fc")
  
  # Set the correct order of columns 
  data_avg_log2fc_genes_of_interest_go_process_long$time_point <- factor(data_avg_log2fc_genes_of_interest_go_process_long$time_point, 
                                                                         levels = c("avg_log2fc_6_hours",
                                                                                    "avg_log2fc_12_hours",
                                                                                    "avg_log2fc_1_days",
                                                                                    "avg_log2fc_3_days",
                                                                                    "avg_log2fc_1_week",
                                                                                    "avg_log2fc_3_weeks"))
  
  # Create the ggplot for glycolysis genes
  go_genes_plot <- ggplot(data_avg_log2fc_genes_of_interest_go_process_long, 
                          aes(x = time_point, y = avg_log2fc, group = gene, color = gene)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_color_brewer(palette = "Set1") +
    xlab("Time Point") +
    ylab("Average log 2 fold change expression") +
    #ggtitle("Average log 2 fold change expression of specific genes") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +
    # theme(axis.title.x = element_blank(), # These three lines excludes the x-axis
    #       axis.text.x = element_blank(),
    #       axis.ticks.x = element_blank()) +
    scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
  
  print(go_genes_plot)
}