# Please help me

# From master script
###########################################################################
# File path to the Excel file
file_path <- "C:/Users/siljeew/snRNAseq/Animal_overview_TP_copy.xlsx"

# Prepare data from original excel file. Renaming columns, excluding some columns, removing space etc. 
processed_data <- process_excel_data(file_path, output_view = TRUE)
View(processed_data) # Check the results

# Define variables of interest and output paths
variable_of_interest <- "hw"  # Heart weight is the desired variable
output_dir <- "C:/Users/siljeew/snRNAseq/Data"
output_csv_path <- file.path(output_dir, "shapiro_wilks_results.csv")

# Ensure the directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Perform Shapiro-Wilk test for normality
perform_shapiro_wilks_test(processed_data, variable_of_interest, output_csv_path)


# Visual investigation of normality: qq-plot
output_dir_distribution <- "Y:/Silje/snRNAseq/Plots/Distribution"

# Run the create_qqplot function

# Function to create and save QQ-plots
create_qqplot <- function(data, variable, output_dir_distribution) {
  # Ensure the output directory exists
  if (!dir.exists(output_dir_distribution)) {
    dir.create(output_dir_distribution, recursive = TRUE)
  }
  
  # Extract unique time points from the dataset
  time_points <- unique(data$time_point)
  
  # Loop through each time point and create qq-plot for each condition
  for (time_point in time_points) {
    # Separate the groups by condition at the current time point
    variable_sham <- data %>%
      filter(condition == 'sham', time_point == time_point) %>%
      select(all_of(variable)) %>%
      unlist() %>%
      na.omit()
    
    variable_ab <- data %>%
      filter(condition == 'AB', time_point == time_point) %>%
      select(all_of(variable)) %>%
      unlist() %>%
      na.omit()
    
    # Debugging prints
    print(paste("Time Point:", time_point))
    print(paste("Sham Data Length:", length(variable_sham)))
    print(paste("AB Data Length:", length(variable_ab)))
    
    # Check the data to ensure proper subsetting
    print(head(variable_sham))
    print(head(variable_ab))
    
    # Create qq-plot for each group
    qq_plot_sham <- ggqqplot(data.frame(variable = variable_sham), x = "variable", ggtheme = theme_bw()) + 
      ggtitle(paste("QQ Plot for sham at time point", time_point)) +
      labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
    
    qq_plot_ab <- ggqqplot(data.frame(variable = variable_ab), x = "variable", ggtheme = theme_bw()) + 
      ggtitle(paste("QQ Plot for AB at time point", time_point)) +
      labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
    
    # Combine the plots
    qq_combined_plot <- ggarrange(qq_plot_sham, qq_plot_ab, ncol = 2, nrow = 1)
    
    # Save the plot
    sanitized_time_point <- gsub("[^a-zA-Z0-9]", "_", time_point)
    file_name <- file.path(output_dir_distribution, paste("qq_plot_", sanitized_time_point, "_", variable, ".png", sep = ""))
    print(file_name)  # Debug print for file name
    
    # Save the plot
    ggsave(file_name, plot = qq_combined_plot, device = "png", width = 40, height = 20, units = "cm")
  }
}



# Run the function
qq_plots <- create_qqplot(processed_data, variable_of_interest, output_dir_distribution)
