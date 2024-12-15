# Boxplot of heart weight
# Input: table (excel?) of animal data
# Output: histograms of all samples, by condition and time_point
# Silje Wilhelmsen

load_libraries <- function() {
  library(tidyverse)
  library(ggpubr)
  library(rstatix)
}

# For testing: this does not make the points appear in the middle of the boxes
create_boxplot <- function(processed_data_highlighted, condition, variable_of_interest) {
  
  # Create plot
  boxplot <- ggplot(processed_data_highlighted, aes(x = !!sym("time_point"), y = !!sym(variable_of_interest)), fill = condition) +
    geom_boxplot(aes(fill = !!sym(condition)), position = position_dodge(width = 0.5), alpha = 0.7) +
    geom_jitter(data = filter(processed_data_highlighted, Used_already_for == "Multiome seq"),
                aes(x = !!sym("time_point"), y = !!sym(variable_of_interest), color = Used_already_for, group = condition),
                size = 2, shape = 21, fill = "red", 
                position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), # Use position_jitterdodge here with small jitter width
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
  
  # Display the plot
  print(boxplot)
}


boxplot <- create_boxplot(processed_data_highlighted, "condition", variable_of_interest)

# ORiginal - does this have centred points? 

create_boxplot <- function(processed_data_highlighted, condition, variable_of_interest) {
  
  # Create plot
  boxplot <- ggplot(processed_data_highlighted, aes(x = !!sym("time_point"), y = !!sym(variable_of_interest)), fill = condition) +
    geom_boxplot(aes(fill = !!sym(condition)), position = position_dodge(width = 0.5), alpha = 0.7) +
    geom_jitter(data = filter(processed_data_highlighted, Used_already_for == "Multiome seq"),
               aes(x = !!sym("time_point"), y = !!sym(variable_of_interest), color = Used_already_for, group = condition),
               size = 2, shape = 21, fill = "red", position = position_dodge(width = 0.5), show.legend = TRUE) +
    
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
  
  # Display the plot
  print(boxplot)
}

boxplot <- create_boxplot(processed_data_highlighted, "condition", variable_of_interest)


# processed_data_highlighted <- processed_data %>%
#   mutate(is_highlighted = ifelse(sample_id %in% highlighted_samples, TRUE, FALSE))
# 
# create_boxplot(processed_data_highlighted, "condition", variable_of_interest)
# 
# 
# ### This is a good one
# create_boxplot <- function(processed_data_highlighted, condition, variable_of_interest) {
#   
#   # Create plot
#   boxplot <- ggplot(processed_data_highlighted, aes_string(x = "time_point", y = variable_of_interest, fill = condition)) +
#     geom_boxplot(position = position_dodge(width = 0.5), alpha = 0.7, show.legend = TRUE) +
#     geom_point(data = filter(processed_data_highlighted, is_highlighted), 
#                aes_string(x = "time_point", y = variable_of_interest, group = condition), 
#                size = 2, shape = 21, fill = "red", 
#                position = position_dodge(width = 0.5), show.legend = FALSE) +
#     labs(title = paste(variable_of_interest, "at all time points"), 
#          x = "Time point", 
#          y = variable_of_interest) +
#     scale_fill_manual(name = "Condition", values = c("SHAM" = "white", "ORAB" = "coral"),
#                       labels = c("SHAM" = "SHAM", "ORAB" = "ORAB")) +
#     theme_minimal() +
#     theme(legend.position = "right") +
#     stat_compare_means(aes(label = ..p.signif..), 
#                        method = "t.test", 
#                        comparisons = list(c("SHAM", "ORAB")), 
#                        label.y = c(1.1, 1.2, 1.3, 1.4, 1.5, 1.6), 
#                        position = position_dodge(width = 0.8))
#   
#   # Display the plot
#   print(boxplot)
# }
# 
# create_boxplot(processed_data_highlighted, "condition", variable_of_interest)
# 



# # Do test and store into something
# 
# p_values <- data.frame(
#   time_point = c("6_hours", "12_hours", "1_day", "3_days", "1_week", "3_weeks"),
#   group1 = "sham", 
#   group2 = "AB",
#   p = c(0.1542, 0.8254, 0.001681, 0.002226, 0.0000448, 0.00001046)
# )
# 
# p_values <- p_values %>%
#   mutate(y.position = c(1.1, 1.2, 1.3, 1.4, 1.5, 1.6))
# 
# View(p_values)
# 
# 
# highlighted_samples <- processed_data %>%
#   filter(Used_already_for == "Multiome seq") %>%
#   pull(sample_id)
# 
# processed_data_highlighted <- processed_data %>%
#   mutate(is_highlighted = ifelse(sample_id %in% highlighted_samples, TRUE, FALSE))
# 
# View(processed_data_highlighted)

# # Create the boxplot with significance annotations
# boxplot_hw <- ggplot(processed_data_highlighted, aes(x = condition, y = hw, fill = condition)) +
#   geom_boxplot(alpha = 0.7, show.legend = FALSE) +
#   geom_point(data = filter(processed_data_highlighted, is_highlighted),
#              aes(x = condition, y = hw),
#              size = 2, shape = 21, fill = "red") +
#   facet_grid(. ~ time_point) +
#   labs(title = "Heart weight (mg) ratio of all samples", 
#        x = "condition", 
#        y = "hw") +
#   scale_fill_manual(name = "condition", values = c("sham" = "white", "AB" = "black"), 
#                     labels = c("sham" = "Sham", "AB" = "AB")) +
#   theme_minimal() +
#   theme(legend.position = "right") +
#   stat_pvalue_manual(p_values, label = "p", tip.length = 0.01)
# 
# # Display the plot
# print(boxplot_hw)
# 
# 
# 
# # Basic boxplot
# simple_boxplot <- ggplot(processed_data_highlighted, aes(x = condition, y = hw, fill = condition)) +
#   geom_boxplot(alpha = 0.7) +
#   facet_grid(.~time_point)
# 
# 
# # Display the simple plot
# print(simple_boxplot)
# 
# 
# #################################
# # boxplot
# boxplot_LVW_BW <- processed_data %>%
#   ggplot(aes(x = condition,  y = `LVW/BW`, fill = condition)) +
#   geom_boxplot(alpha = 0.7) +
#   facet_grid(. ~ time_point) +
#   labs(title = "Left ventricular eight/body weight ratio of all samples",
#        x = "LVW/BW", 
#        y = "Frequency") +
#   scale_fill_manual(values = c("sham" = "white", "AB" = "black")) + 
#   theme_minimal() +
#   theme(legend.position = "none")
# 
# print(boxplot_LVW_BW)
# 
# # Save
# 
# 
# # With specific highlighting
# animal_data_rnaseq_highlighted <- animal_data_rnaseq_renamed %>%
#   mutate(is_highlighted = ifelse(sample_id %in% highlighted_samples, TRUE, FALSE))
# 
# # Create a boxplot with highlighted points
# boxplot_LVW_BW <- animal_data_rnaseq_highlighted %>%
#   ggplot(aes(x = condition, y = `LVW/BW`, fill = condition)) +
#   geom_boxplot(alpha = 0.7, show.legend = FALSE) +
#   geom_point(data = filter(animal_data_rnaseq_highlighted, is_highlighted),
#              aes(x = condition, y = `LVW/BW`, color = "Samples for RNA sequencing"),
#              size = 2, shape = 21, fill = "red", show.legend = TRUE) +
#   facet_grid(. ~ time_point) +
#   labs(title = "Left ventricular eight/body weight ratio of all samples", 
#        x = "Condition", 
#        y = "LVW/BW") +
#   scale_fill_manual(values = c("sham" = "white", "AB" = "black"), guide = "none") + 
#   theme_minimal() +
#   theme(legend.position = "right")
# 
# # Display the plot
# print(boxplot_LVW_BW)
# 
# 
# 
