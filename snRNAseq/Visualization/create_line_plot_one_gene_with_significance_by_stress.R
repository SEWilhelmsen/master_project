# Create line plot with significance labels
# Silje Wilhelmsen



# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary$Timepoint <- factor(go_process_summary$Timepoint, 
                                       levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary$Stress_Status <- factor(go_process_summary$Stress_Status,
                                           levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))

# Create plot for BDH1
#################################################################################################

go_process_summary_bdh1 <- go_process_summary

# Filter just for significant labels; assume y_position is aligned
significant_labels <- go_process_summary_bdh1 %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = percentage_to_sham - 1) # Adjust y_position as needed


# Filter just for significant labels; assume y_position is aligned
stressed_not_stressed <- go_process_summary_bdh1 %>%
  filter(stress_vs_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 15) # Adjust y_position as needed


# Include ANOVA value? 

single_gene_plot_bdh1 <- ggplot(go_process_summary_bdh1, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(size = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression compared to SHAM CM (%)", 
       title = paste("BDH1"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) + # Check that factor are set
  scale_y_continuous(breaks = seq(0, 300, by = 40), limits = c(60, 300)) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 18, color = "black", margin = margin(t = 15)), 
        axis.text.y = element_text(size = 16, colour = "black"), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 26, margin = margin(r = 15)), # Change margins: r = right, b = bottom, t = top, l = left
        plot.title = element_text(size = 26, margin = margin(b = 20)),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.line.y.right = element_blank()) +
  
  geom_text(data = significant_labels, aes(x = Timepoint, y = y_position, label = significance_label), 
            size = 8, vjust = -0.5, color = "black") +
  geom_text(data = stressed_not_stressed, aes(x = Timepoint, y = y_position, label = stress_vs_not_stressed), 
            size = 6, vjust = -0.5, color = "black")

single_gene_plot_bdh1




# Create plot for OXCT1
################################################################################################
# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary$Timepoint <- factor(go_process_summary$Timepoint, 
                                       levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary$Stress_Status <- factor(go_process_summary$Stress_Status,
                                           levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))


# # Filter just for significant labels; assume y_position is aligned
# significant_labels <- go_process_summary %>%
#   filter(significance_label != "ns") %>%
#   mutate(y_position = percentage_to_sham + 0.3) # Adjust y_position as needed

go_process_summary_oxct1 <- go_process_summary

# Modify specific timepoints if necessary
significant_labels <- go_process_summary_oxct1 %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = ifelse(Timepoint %in% c("12 Hours", "3 Weeks") & Stress_Status %in% "Stressed CM", 
                             percentage_to_sham + 1, ## Custom adjustment for overlapping points
                             percentage_to_sham - 0.7))  # General adjustment


# Filter just for significant labels; assume y_position is aligned
stressed_not_stressed <- go_process_summary_oxct1 %>%
  filter(stress_vs_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 2) # Adjust y_position as needed

# # Modify specific timepoints if necessary
# stressed_not_stressed <- go_process_summary %>%
#   filter(stress_vs_not_stressed != "ns") %>%
#   mutate(y_position = ifelse(Timepoint %in% c("SpecificTimepoint1", "SpecificTimepoint2"), 
#                              percentage_to_sham + 4.6, ## Custom adjustment for overlapping points
#                              percentage_to_sham + 4))  # General adjustment

single_gene_plot_oxct1 <- ggplot(go_process_summary_oxct1, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression compared to SHAM CM (%)", 
       title = paste("OXCT1"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 160, by = 10), limits = c(70, 120)) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 18, color = "black", margin = margin(t = 15)), 
        axis.text.y = element_text(size = 16, colour = "black"), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 26, margin = margin(r = 15)), # Change margins: r = right, b = bottom, t = top, l = left
        plot.title = element_text(size = 26, margin = margin(b = 20)),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.line.y.right = element_blank()) +
  
  geom_text(data = significant_labels, aes(x = Timepoint, y = y_position, label = significance_label), 
            size = 8, vjust = -0.5, color = "black") +
  geom_text(data = stressed_not_stressed, aes(x = Timepoint, y = y_position, label = stress_vs_not_stressed), 
            size = 6, vjust = -0.5, color = "black")

single_gene_plot_oxct1








# Make a figure consisting of multiple figures
##############################################################################################

go_plot_hydroxybutyrate_dehydrogenase <- go_plot_hydroxybutyrate_dehydrogenase +
  theme(plot.margin = unit(c(1,1,2,2.5), "cm"))  # top, right, bottom, left margin
go_plot_ketone_catabolism <- go_plot_ketone_catabolism +
  theme(plot.margin = unit(c(1,1,2,2.5), "cm"))  # top, right, bottom, left margin

gene_plot_bdh1 <- gene_plot_bdh1 +
  theme(plot.margin = unit(c(1,1,2,2.5), "cm"))  # top, right, bottom, left margin
gene_plot_oxct1_acat <- gene_plot_oxct1_acat +
  theme(plot.margin = unit(c(1,1,2,2.5), "cm"))  # top, right, bottom, left margin

single_gene_plot_bdh1 <- single_gene_plot_bdh1 +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin
single_gene_plot_oxct1 <- single_gene_plot_oxct1 +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin


# Example of laying out in a 3x3 grid
total_ketone_plot <- (go_plot_hydroxybutyrate_dehydrogenase | go_plot_ketone_catabolism) /
  (gene_plot_bdh1 | gene_plot_oxct1_acat) /
  (single_gene_plot_bdh1 | single_gene_plot_oxct1)

# Print the combined plot
print(total_ketone_plot)

# Save the combined plot—adjust width and height for space between plots
ggsave(file.path(output_dir_plot, "total_ketone_plot.pdf"), plot = total_ketone_plot, width = 24, height = 28)

ggsave(file.path(output_dir_plot, "total_ketone_plot.png"), plot = total_ketone_plot, width = 24, height = 28)


