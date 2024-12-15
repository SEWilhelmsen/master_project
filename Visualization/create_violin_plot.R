# Create violin plot
# Silje Wilhelmsen

load_libraries <- function() {
  library(ggplot2)
}

# Function to create violin plot
create_violin_plot <- function(combined_data_long, condition, gene) {
  
  # Set the correct order 
  combined_data_long$time_point <- factor(combined_data_long$time_point,
                                    levels = c("3 Days",
                                               "1 week",
                                               "3 Weeks"))
  
  # Create plot
  violin_plot <- ggplot(combined_data_long, aes(x = time_point, y = expression, fill = condition)) +
    geom_violin(position = position_dodge(width = 0.5), alpha = 0.7, show.legend = TRUE) +
    labs(title = paste("Expression of ", gene, " over time"), 
         x = "Time point", 
         y = paste("Expression of ", gene)) +
    scale_fill_manual(name = "condition", values = c("SHAM" = "white", "AB" = "coral"), 
                      labels = c("sham" = "SHAM", "AB" = "ORAB")) +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Display the plot
  print(violin_plot)
}