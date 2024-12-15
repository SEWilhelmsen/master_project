### QQ-plot to investigate distribution of a variable (ikke ferdig)
# Input: processed_data
# Output: qq-plot of variable_of_interest for all time points, sham vs AB
## Silje Wilhelmsen

# # Load libraries
# library(dplyr)
# library(ggpubr)
# library(ggplot2)


# Function to create and save QQ-plots
create_qqplot <- function(processed_data, variable, output_dir_distribution) {
  # Ensure the output directory exists
  if (!dir.exists(output_dir_distribution)) {
    dir.create(output_dir_distribution, recursive = TRUE)
  }
  
  # Extract unique time points from the dataset
  time_points <- unique(processed_data$time_point)
  
  # Loop through each time point and create qq-plot for each condition
  for (time_point in time_points) {
    # Separate the groups by condition at the current time point
    variable_sham <- processed_data %>%
      filter(condition == 'sham', time_point == time_point) %>%
      select(all_of(variable)) %>%
      unlist() %>%
      na.omit()
    
    variable_ab <- processed_data %>%
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




# Function to create and save histogram
create_histogram <- function(processed_data, variable_of_interest, output_dir, binwidth = 0.3, fill_values = c("sham" = "white", "AB" = "black")) {
  
  # Create the histogram plot
  histogram_plot <- processed_data %>%
    ggplot(aes_string(x = variable_of_interest, fill = "condition")) +
    geom_histogram(binwidth = binwidth, color = "black", alpha = 0.7, position = "dodge") +
    facet_grid(condition ~ time_point) +
    labs(title = paste("Histogram of", variable_of_interest, "by Condition and Time Point"), 
         x = variable_of_interest, 
         y = "Frequency") +
    scale_fill_manual(values = fill_values) + 
    theme_minimal() +
    theme(legend.position = "none")
  
  # Print the plot
  print(histogram_plot)
  
  # Save the plot as a PNG file
  output_file_path <- file.path(output_dir, paste("histogram_", variable_of_interest, ".png"))
  ggsave(output_file_path, plot = histogram_plot, width = 10, height = 6)
}