# Create bar plot of one gene ontology process
# Silje Wilhelmsen



# Function to generate plots with significance stars
generate_plot <- function(summary_data, significance_data, go_process_of_interest) {
  View(significance_data)
  
  significance_data <- significance_data %>%
    filter(process == go_process_of_interest) %>%
    dplyr::select(time_point, condition, p_adj, significance_label) %>%
    filter(significance_label != "ns")
  
  
  significance_data <- significance_data %>%
    left_join(summary_data %>%
                group_by(time_point) %>%
                summarize(max_mean_expression = max(mean_expression, na.rm = TRUE)),
              by = "time_point") %>%
    mutate(y.position = max_mean_expression * 1.01) %>%
    dplyr::select(time_point, condition, p_adj, significance_label, y.position)
  
  print(significance_data)
  
  summary_data$time_point <- factor(summary_data$time_point,
                                    levels = c("6 Hours",
                                               "12 Hours",
                                               "1 Day",
                                               "3 Days",
                                               "1 Week",
                                               "3 Weeks"))
  
  max_y_position <- max(significance_data$y.position, na.rm = TRUE)
  
  go_plot <- ggplot(summary_data, aes(x = time_point, y = mean_expression, fill = condition)) +
    geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
    labs(x = NULL, 
         y = "Mean Expression Level",
         title = "Fatty Acid Transport",
         fill = "Condition") +
    theme_classic() +
    scale_fill_manual(values = c("SHAM" = "grey22", "ORAB" = "coral")) +
    scale_y_continuous(limits = c(0, max_y_position + 0.08)) +  # Add extra space above the bars
    theme(axis.text.x = element_text(hjust = 0.5, size = 18, color = "black", margin = margin(t = 15)), 
          axis.text.y = element_text(size = 26, colour = "black"), 
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 26, margin = margin(r = 15)), # Change margins: r = right, b = bottom, t = top, l = left
          plot.title = element_text(size = 26, margin = margin(b = 20)),
          legend.title = element_text(size = 30),
          legend.text = element_text(size = 30),
          axis.line.y.right = element_blank()) +
    geom_text(data = significance_data,
              aes(x = time_point, y = y.position, label = significance_label, fill = NULL),
              size = 8, vjust = -0.5, color = "black", position = position_dodge(width = 0.75)) 
  
  
  return(go_plot)
}
 
go_plot_fatty_acid_uptake <- generate_plot(go_process_summary, go_process_summary, go_process_of_interest)
print(go_plot_fatty_acid_uptake)
go_process_summary_fatty_acid_uptake <- go_process_summary

# go_plot_b_oxidation <- generate_plot(go_process_summary, go_process_summary, go_process_of_interest)
# print(go_plot_b_oxidation)
# go_process_summary_b_oxidation <- go_process_summary





# ggsave(file.path(output_dir_plot, paste(go_process_of_interest, ".png", sep = "")), plot = go_plot_cardiac_muscle_growth, width = 8, height = 6)
# ggsave(file.path(output_dir_plot, paste(go_process_of_interest, ".pdf", sep = "")), plot = go_plot_cardiac_muscle_growth, width = 8, height = 6)

