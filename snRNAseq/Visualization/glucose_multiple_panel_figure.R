# Create line plot with significance labels
# Silje Wilhelmsen




# Create plot for slc2a4
#################################################################################################
# Load data
# go_process_summary_slc2a4 <- read_excel("Data/Stress_Status/go_process_summary_slc2a4.xlsx")
go_process_summary_slc2a4 <- read_excel("C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/go_process_summary_t_test_slc2a4.xlsx")

# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary_slc2a4$Timepoint <- factor(go_process_summary_slc2a4$Timepoint, 
                                       levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary_slc2a4$Stress_Status <- factor(go_process_summary_slc2a4$Stress_Status,
                                           levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))


# Modify specific timepoints if necessary
significant_labels <- go_process_summary_slc2a4 %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = ifelse(Timepoint %in% c("6 Hours", "12 Hours") & Stress_Status %in% c("Not stressed CM"),
                             percentage_to_sham - 0, 
                             percentage_to_sham - 5))  # General adjustment


# Modify specific timepoints if necessary
stressed_not_stressed <- go_process_summary_slc2a4 %>%
  filter(stressed_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 8)  # General adjustment


single_gene_plot_slc2a4 <- ggplot(go_process_summary_slc2a4, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression (%)", 
       title = paste("SLC2A4 Expression Compared to SHAM"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 150, by = 20), limits = c(30, 150)) +
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

single_gene_plot_slc2a4




# # Create plot for slc27a1
# ################################################################################################
# # Load data
# # go_process_summary_slc27a1 <- read_excel("Data/Stress_Status/go_process_summary_slc27a1.xlsx")
# go_process_summary_slc27a1 <- read_excel("C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/go_process_summary_t_test_slc27a1.xlsx")
# 
# 
# # Setting the Timepoint as a factor to ensure order in the plot
# go_process_summary_slc27a1$Timepoint <- factor(go_process_summary_slc27a1$Timepoint, 
#                                               levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
# 
# # Order the Stress_Status levels
# go_process_summary_slc27a1$Stress_Status <- factor(go_process_summary_slc27a1$Stress_Status,
#                                                   levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))
# 
# 
# # Modify specific timepoints if necessary
# significant_labels <- go_process_summary_slc27a1 %>%
#   filter(significance_label != "ns") %>%
#   mutate(y_position = ifelse(Timepoint %in% c("3 Weeks") & Stress_Status %in% c("Not stressed CM"), 
#                              percentage_to_sham - 4, 
#                              percentage_to_sham - 5))  # General adjustment
# 
# 
# # Filter just for significant labels; assume y_position is aligned
# stressed_not_stressed <- go_process_summary_slc27a1 %>%
#   filter(stressed_not_stressed != "ns") %>%
#   mutate(y_position = ifelse(Timepoint %in% c("12 Hours"), 
#                              percentage_to_sham + 10, 
#                              percentage_to_sham + 2)) # Adjust y_position as needed
# 
# 
# single_gene_plot_slc27a1 <- ggplot(go_process_summary_slc27a1, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
#   geom_line(linewidth = 1.5) +  
#   geom_point(size = 3) +  
#   labs(x = NULL, 
#        y = "Expression (%)", 
#        title = paste("SLC27A1 Expression Compared to SHAM"),
#        color = "Stress status") +  # Grouping title
#   theme_classic() +
#   scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
#   scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
#   scale_y_continuous(breaks = seq(0, 110, by = 20), limits = c(40, 110)) +
#   theme(axis.text.x = element_text(hjust = 0.5, size = 18, color = "black", margin = margin(t = 15)), 
#         axis.text.y = element_text(size = 20, colour = "black"), 
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size = 26, margin = margin(r = 15)), # Change margins: r = right, b = bottom, t = top, l = left
#         plot.title = element_text(size = 26, margin = margin(b = 20)),
#         legend.title = element_text(size = 30),
#         legend.text = element_text(size = 30),
#         axis.line.y.right = element_blank()) +
#   
#   geom_text(data = significant_labels, aes(x = Timepoint, y = y_position, label = significance_label), 
#             size = 8, vjust = -0.5, color = "black") +
#   geom_text(data = stressed_not_stressed, aes(x = Timepoint, y = y_position, label = stressed_not_stressed), 
#             size = 6, vjust = -0.5, color = "black")
# 
# single_gene_plot_slc27a1

# ggsave(file.path(output_dir_plot, paste("slc27a1.pdf", sep = "")), plot = single_gene_plot_slc27a1, width = 8, height = 6)



# Create plot for GCK 
################################################################################################
# Load data
# go_process_summary_gck <- read_excel("Data/Stress_Status/go_process_summary_gck.xlsx")
go_process_summary_gck <- read_excel("C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/go_process_summary_t_test_gck.xlsx")


# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary_gck$Timepoint <- factor(go_process_summary_gck$Timepoint, 
                                               levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary_gck$Stress_Status <- factor(go_process_summary_gck$Stress_Status,
                                                   levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))

# Modify specific timepoints if necessary
significant_labels <- go_process_summary_gck %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = ifelse(Timepoint %in% c("12 Hours") & Stress_Status %in% c("Not stressed CM"), 
                             percentage_to_sham - 0,
                             percentage_to_sham - 12))  # General adjustment


# Filter just for significant labels; assume y_position is aligned
stressed_not_stressed <- go_process_summary_gck %>%
  filter(stressed_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 50) # Adjust y_position as needed


single_gene_plot_gck <- ggplot(go_process_summary_gck, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression (%)", 
       title = paste("GCK Expression Compared to SHAM"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 260, by = 30), limits = c(20, 260)) +
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

single_gene_plot_gck



# Create plot for HK1
################################################################################################
# Load data
# go_process_summary_hk1 <- read_excel("Data/Stress_Status/go_process_summary_hk1.xlsx")
go_process_summary_hk1 <- read_excel("C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/go_process_summary_t_test_hk1.xlsx")


# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary_hk1$Timepoint <- factor(go_process_summary_hk1$Timepoint, 
                                               levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary_hk1$Stress_Status <- factor(go_process_summary_hk1$Stress_Status,
                                                   levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))


# Modify specific timepoints if necessary
significant_labels <- go_process_summary_hk1 %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = ifelse(Timepoint %in% c("1 Day"),
                             percentage_to_sham - 8, 
                             percentage_to_sham - 5))  # General adjustment


# Filter just for significant labels; assume y_position is aligned
stressed_not_stressed <- go_process_summary_hk1 %>%
  filter(stressed_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 15) # Adjust y_position as needed


single_gene_plot_hk1 <- ggplot(go_process_summary_hk1, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression (%)", 
       title = paste("HK1 Expression Compared to SHAM"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 220, by = 20), limits = c(70, 220)) +
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

single_gene_plot_hk1



# Create plot for HK2
################################################################################################
# Load data
# go_process_summary_hk2 <- read_excel("Data/Stress_Status/go_process_summary_hk2.xlsx")
go_process_summary_hk2 <- read_excel("C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/go_process_summary_t_test_hk2.xlsx")


# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary_hk2$Timepoint <- factor(go_process_summary_hk2$Timepoint, 
                                           levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary_hk2$Stress_Status <- factor(go_process_summary_hk2$Stress_Status,
                                               levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))



# Modify specific timepoints if necessary
significant_labels <- go_process_summary_hk2 %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = percentage_to_sham - 1.5)  # General adjustment


# Filter just for significant labels; assume y_position is aligned
stressed_not_stressed <- go_process_summary_hk2 %>%
  filter(stressed_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 6) # Adjust y_position as needed


single_gene_plot_hk2 <- ggplot(go_process_summary_hk2, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression (%)", 
       title = paste("HK2 Expression Compared to SHAM"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 150, by = 20), limits = c(80, 150)) +
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

single_gene_plot_hk2



# Create plot for PFKM
################################################################################################
# Load data
# go_process_summary_pfkm <- read_excel("Data/Stress_Status/go_process_summary_pfkm.xlsx")
go_process_summary_pfkm <- read_excel("C:/Users/siljeew/Master_project/snRNAseq/Data/Stress_Status/go_process_summary_t_test_pfkm.xlsx")


# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary_pfkm$Timepoint <- factor(go_process_summary_pfkm$Timepoint, 
                                           levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary_pfkm$Stress_Status <- factor(go_process_summary_pfkm$Stress_Status,
                                               levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))


# Modify specific timepoints if necessary
significant_labels <- go_process_summary_pfkm %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = percentage_to_sham - 3)  # General adjustment

# Filter just for significant labels; assume y_position is aligned
stressed_not_stressed <- go_process_summary_pfkm %>%
  filter(stressed_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 3) # Adjust y_position as needed


single_gene_plot_pfkm <- ggplot(go_process_summary_pfkm, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression (%)", 
       title = paste("PFKM Expression Compared to SHAM"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 140, by = 20), limits = c(50, 140)) +
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

single_gene_plot_pfkm




# To combine plots:
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/glucose_multiple_panel_figure_combined.R")

