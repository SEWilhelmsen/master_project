# Create line plot with significance labels
# Silje Wilhelmsen

# Necessary preparation:
# file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/run_create_barplot_go_process_by_stress.R") #This is not finished
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/run_create_gene_plot_by_stress_.R")



# Create plot for IDH1
################################################################################################
# Load data
# go_process_summary_idh1 <- read_excel("Data/Stress_Status/go_process_summary_idh1.xlsx")
go_process_summary_idh1 <- read_excel("Data/Stress_Status/go_process_summary_t_test_idh1.xlsx")


# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary_idh1$Timepoint <- factor(go_process_summary_idh1$Timepoint, 
                                       levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary_idh1$Stress_Status <- factor(go_process_summary_idh1$Stress_Status,
                                           levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))


# Modify specific timepoints if necessary
significant_labels <- go_process_summary_idh1 %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = percentage_to_sham - 2.5)  # General adjustment


# Filter just for significant labels; assume y_position is aligned
stressed_not_stressed <- go_process_summary_idh1 %>%
  filter(stressed_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 2.5) # Adjust y_position as needed


single_gene_plot_idh1 <- ggplot(go_process_summary_idh1, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression (%)", 
       title = paste("IDH1 Expression Compared to SHAM"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 140, by = 20), limits = c(40, 140)) +
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

single_gene_plot_idh1



# Create plot for IDH2
################################################################################################
# Load data
# go_process_summary_idh2 <- read_excel("Data/Stress_Status/go_process_summary_idh2.xlsx")
go_process_summary_idh2 <- read_excel("Data/Stress_Status/go_process_summary_t_test_idh2.xlsx")


# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary_idh2$Timepoint <- factor(go_process_summary_idh2$Timepoint, 
                                       levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary_idh2$Stress_Status <- factor(go_process_summary_idh2$Stress_Status,
                                           levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))

#go_process_summary_idh2 <- go_process_summary

# Modify specific timepoints if necessary
significant_labels <- go_process_summary_idh2 %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = percentage_to_sham - 1)  # General adjustment


# Filter just for significant labels; assume y_position is aligned
stressed_not_stressed <- go_process_summary_idh2 %>%
  filter(stressed_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 4) # Adjust y_position as needed


single_gene_plot_idh2 <- ggplot(go_process_summary_idh2, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression (%)", 
       title = paste("IDH2 Expression Compared to SHAM"),
       color = "Stress status") +  
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 110, by = 20), limits = c(40, 110)) +
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

single_gene_plot_idh2



# Create plot for IDH3B
################################################################################################
# Load data
# go_process_summary_idh3b <- read_excel("Data/Stress_Status/go_process_summary_idh3b.xlsx")
go_process_summary_idh3b <- read_excel("Data/Stress_Status/go_process_summary_t_test_idh3b.xlsx")


# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary_idh3b$Timepoint <- factor(go_process_summary_idh3b$Timepoint, 
                                       levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary_idh3b$Stress_Status <- factor(go_process_summary_idh3b$Stress_Status,
                                           levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))


# Modify specific timepoints if necessary
significant_labels <- go_process_summary_idh3b %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = percentage_to_sham - 2)  


# Filter just for significant labels; assume y_position is aligned
stressed_not_stressed <- go_process_summary_idh3b %>%
  filter(stressed_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 4) 


single_gene_plot_idh3b <- ggplot(go_process_summary_idh3b, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression (%)", 
       title = paste("IDH3B Expression Compared to SHAM"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 160, by = 20), limits = c(50, 160)) +
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

single_gene_plot_idh3b



# To combine plots:
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/tca_multiple_panel_figure_combined.R")
