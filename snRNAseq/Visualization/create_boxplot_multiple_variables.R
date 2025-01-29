# Create multiple boxplots in a combined plot
# Silje Wilhelmsen

load_libraries <- function() {
  library(tidyverse)
  library(ggplot2)
  library(readxl)
  library(ggpubr)
  library(rstatix)
  library(gridExtra)
}


# Define the function to create a boxplot
create_boxplot <- function(processed_data_highlighted, condition, variable_of_interest) {
  
  # Create plot
  boxplot <- ggplot(processed_data_highlighted, aes(x = !!sym("time_point"), y = !!sym(variable_of_interest)), fill = condition) +
    geom_boxplot(aes(fill = !!sym(condition)), position = position_dodge(width = 0.5), alpha = 0.7) +
    geom_jitter(data = filter(processed_data_highlighted, Used_already_for == "Multiome seq"),
                aes(x = !!sym("time_point"), y = !!sym(variable_of_interest), color = Used_already_for, group = condition),
                size = 2, shape = 21, fill = "red", 
                position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
                show.legend = TRUE) +
    
    labs(title = paste(variable_of_interest, "at all time points"), 
         x = "Time point", 
         y = variable_of_interest) +
    
    scale_fill_manual(name = "Condition", values = c("SHAM" = "white", "ORAB" = "coral"),
                      labels = c("SHAM" = "SHAM", "ORAB" = "ORAB")) +
    
    scale_color_manual(name = "Sample", values = c("Multiome seq" = "red"), labels = c("RNA sequencing sample")) +
    
    theme_minimal() +
    theme(legend.position = "right") +
    guides(fill = guide_legend(override.aes = list(shape = NA)),
           color = guide_legend(override.aes = list(shape = 21, fill = "red")))
  
  return(boxplot)
}

