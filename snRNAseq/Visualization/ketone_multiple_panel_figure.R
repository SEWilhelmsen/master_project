# Create line plot with significance labels
# Silje Wilhelmsen

# Necessary preparation:
# file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/run_create_barplot_go_process_by_stress.R") #This is not finished
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/run_create_gene_plot_by_stress_.R")



# Create plot for BDH1
################################################################################################
# Load data
# go_process_summary_bdh1 <- read_excel("Data/Stress_Status/go_process_summary_bdh1.xlsx")
go_process_summary_bdh1 <- read_excel("Data/Stress_Status/go_process_summary_t_test_bdh1.xlsx")


# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary_bdh1$Timepoint <- factor(go_process_summary_bdh1$Timepoint, 
                                            levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary_bdh1$Stress_Status <- factor(go_process_summary_bdh1$Stress_Status,
                                                levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))

#go_process_summary_bdh1 <- go_process_summary

# Modify specific timepoints if necessary
significant_labels <- go_process_summary_bdh1 %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = percentage_to_sham - 2.5)  # General adjustment


# Filter just for significant labels; assume y_position is aligned
stressed_not_stressed <- go_process_summary_bdh1 %>%
  filter(stressed_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 10) # Adjust y_position as needed


single_gene_plot_bdh1 <- ggplot(go_process_summary_bdh1, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression (%)", 
       title = paste("BDH1 Expression Compared to SHAM"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 280, by = 20), limits = c(60, 280)) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 18, color = "black", margin = margin(t = 15)), 
        axis.text.y = element_text(size = 20, colour = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 26, margin = margin(r = 15)), # Change margins: r = right, b = bottom, t = top, l = left
        plot.title = element_text(size = 26, margin = margin(b = 20)),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30),
        axis.line.y.right = element_blank()) +
  geom_text(data = significant_labels, aes(x = Timepoint, y = y_position, label = significance_label), 
            size = 8, vjust = -0.5, color = "black") +
  geom_text(data = stressed_not_stressed, aes(x = Timepoint, y = y_position, label = stressed_not_stressed), 
            size = 6, vjust = -0.5, color = "black")

single_gene_plot_bdh1



# Create plot for OXCT1
################################################################################################
# Load data
# go_process_summary_oxct1 <- read_excel("Data/Stress_Status/go_process_summary_oxct1.xlsx")
go_process_summary_oxct1 <- read_excel("Data/Stress_Status/go_process_summary_t_test_oxct1.xlsx")


# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary_oxct1$Timepoint <- factor(go_process_summary_oxct1$Timepoint, 
                                            levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary_oxct1$Stress_Status <- factor(go_process_summary_oxct1$Stress_Status,
                                                levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))


# Modify specific timepoints if necessary
significant_labels <- go_process_summary_oxct1 %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = ifelse(Timepoint %in% c("12 Hours", "3 Weeks") & Stress_Status %in% c("Not stressed CM"), 
                             percentage_to_sham + 1,
                             percentage_to_sham - 1)) 


# Filter just for significant labels; assume y_position is aligned
stressed_not_stressed <- go_process_summary_oxct1 %>%
  filter(stressed_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 3)


single_gene_plot_oxct1 <- ggplot(go_process_summary_oxct1, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression (%)", 
       title = paste("OXCT1 Expression Compared to SHAM"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 120, by = 20), limits = c(60, 120)) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 18, color = "black", margin = margin(t = 15)), 
        axis.text.y = element_text(size = 20, colour = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 26, margin = margin(r = 15)), # Change margins: r = right, b = bottom, t = top, l = left
        plot.title = element_text(size = 26, margin = margin(b = 20)),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30),
        axis.line.y.right = element_blank()) +
  
  geom_text(data = significant_labels, aes(x = Timepoint, y = y_position, label = significance_label), 
            size = 8, vjust = -0.5, color = "black") +
  geom_text(data = stressed_not_stressed, aes(x = Timepoint, y = y_position, label = stressed_not_stressed), 
            size = 6, vjust = -0.5, color = "black")

single_gene_plot_oxct1




# To combine plots:
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/ketone_multiple_panel_figure_combined.R")
