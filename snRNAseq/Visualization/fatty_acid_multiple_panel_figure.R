# Create line plot with significance labels
# Silje Wilhelmsen

# # For significance label: 
# t_test_gene_by_stress <- read_excel("Data/Stress_Status/t_test_gene_by_stress.xlsx")
# View(t_test_gene_by_stress)
# 
# # If necessary:
# t_test_gene_by_stress <- t_test_gene_by_stress %>%
#   mutate(significance_label = ifelse(p_adj < 0.001, "***",
#                                      ifelse(p_adj < 0.01, "**",
#                                             ifelse(p_adj < 0.05, "*", "ns")))) %>%
#   mutate(stressed_not_stressed = ifelse(group1 != "SHAM CM" & group2 != "SHAM CM" & p_adj < 0.05, "\u25B2", "ns"))



# Create plot for CD36
#################################################################################################
# Load data
# go_process_summary_cd36 <- read_excel("Data/Stress_Status/go_process_summary_cd36.xlsx")
go_process_summary_cd36 <- read_excel("C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/go_process_summary_t_test_cd36.xlsx")

# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary_cd36$Timepoint <- factor(go_process_summary_cd36$Timepoint, 
                                       levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary_cd36$Stress_Status <- factor(go_process_summary_cd36$Stress_Status,
                                           levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))



# Modify specific timepoints if necessary
significant_labels <- go_process_summary_cd36 %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = ifelse(Timepoint %in% c("1 Week") & Stress_Status %in% "Stressed CM", 
                             percentage_to_sham + 1, ## Custom adjustment for overlapping points
                             percentage_to_sham - 2))  # General adjustment



# Modify specific timepoints if necessary
stressed_not_stressed <- go_process_summary_cd36 %>%
  filter(stressed_not_stressed != "ns") %>%
  mutate(y_position = ifelse(Timepoint %in% c("3 Days", "3 Weeks"), percentage_to_sham + 4, 
                             ifelse(Timepoint %in% c("1 Day"), percentage_to_sham + 3, ## Custom adjustment for overlapping points
                                    percentage_to_sham + 8)))  # General


single_gene_plot_cd36 <- ggplot(go_process_summary_cd36, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression (%)", 
       title = paste("CD36 Expression Compared to SHAM"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 130, by = 10), limits = c(70, 130)) +
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

single_gene_plot_cd36




# Create plot for CPT1A
################################################################################################
# Load data
# go_process_summary_cpt1a <- read_excel("Data/Stress_Status/go_process_summary_cpt1a.xlsx")
go_process_summary_cpt1a <- read_excel("C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/go_process_summary_t_test_cpt1a.xlsx")

# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary_cpt1a$Timepoint <- factor(go_process_summary_cpt1a$Timepoint, 
                                       levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary_cpt1a$Stress_Status <- factor(go_process_summary_cpt1a$Stress_Status,
                                           levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))


# Modify specific timepoints if necessary
significant_labels <- go_process_summary_cpt1a %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = percentage_to_sham - 3)  # General adjustment


# Filter just for significant labels; assume y_position is aligned
stressed_not_stressed <- go_process_summary_cpt1a %>%
  filter(stressed_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 2) # Adjust y_position as needed


single_gene_plot_cpt1a <- ggplot(go_process_summary_cpt1a, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression (%)", 
       title = paste("CPT1A Expression Compared to SHAM"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  # scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 160, by = 20), limits = c(20, 160)) +
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

single_gene_plot_cpt1a



# Create plot for CPT1B
################################################################################################
# Load data
# go_process_summary_cpt1b <- read_excel("Data/Stress_Status/go_process_summary_cpt1b.xlsx")
go_process_summary_cpt1b <- read_excel("C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/go_process_summary_t_test_cpt1b.xlsx")


# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary_cpt1b$Timepoint <- factor(go_process_summary_cpt1b$Timepoint, 
                                       levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary_cpt1b$Stress_Status <- factor(go_process_summary_cpt1b$Stress_Status,
                                           levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))


# Modify specific timepoints if necessary
significant_labels <- go_process_summary_cpt1b %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = percentage_to_sham - 1.5)  # General adjustment


# Filter just for significant labels; assume y_position is aligned
stressed_not_stressed <- go_process_summary_cpt1b %>%
  filter(stressed_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 4) # Adjust y_position as needed


single_gene_plot_cpt1b <- ggplot(go_process_summary_cpt1b, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression (%)", 
       title = paste("CPT1B Expression Compared to SHAM"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 120, by = 10), limits = c(60, 120)) +
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

single_gene_plot_cpt1b



# Create plot for CPT2
################################################################################################
# Load data
# go_process_summary_cpt2 <- read_excel("Data/Stress_Status/go_process_summary_cpt2.xlsx")
go_process_summary_cpt2 <- read_excel("C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/go_process_summary_t_test_cpt2.xlsx")


# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary_cpt2$Timepoint <- factor(go_process_summary_cpt2$Timepoint, 
                                       levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary_cpt2$Stress_Status <- factor(go_process_summary_cpt2$Stress_Status,
                                           levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))



# Modify specific timepoints if necessary
significant_labels <- go_process_summary_cpt2 %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = percentage_to_sham - 3)  # General adjustment


# Filter just for significant labels; assume y_position is aligned
stressed_not_stressed <- go_process_summary_cpt2 %>%
  filter(stressed_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 9) # Adjust y_position as needed


single_gene_plot_cpt2 <- ggplot(go_process_summary_cpt2, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression (%)", 
       title = paste("CPT2 Expression Compared to SHAM"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 160, by = 20), limits = c(20, 160)) +
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

single_gene_plot_cpt2


# Combine plots:
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/fatty_acid_multiple_panel_figure_combined.R")
