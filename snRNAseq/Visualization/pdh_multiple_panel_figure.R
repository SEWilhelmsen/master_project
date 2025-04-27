# Create line plot with significance labels
# Silje Wilhelmsen


# For preparing t-test labels
# file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/prepare_t_test_significance_label_gene_by_stress.R")

# Create plot for PDHA1
#################################################################################################
# Load data
# go_process_summary_pdha1 <- read_excel("Data/Stress_Status/go_process_summary_pdha1.xlsx")
go_process_summary_pdha1 <- read_excel("Data/Stress_Status/go_process_summary_t_test_pdha1.xlsx")



# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary_pdha1$Timepoint <- factor(go_process_summary_pdha1$Timepoint, 
                                       levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary_pdha1$Stress_Status <- factor(go_process_summary_pdha1$Stress_Status,
                                           levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))

# Modify specific timepoints if necessary
significant_labels <- go_process_summary_pdha1 %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = percentage_to_sham - 2)  # General adjustment

# Filter just for significant labels; assume y_position is aligned
stressed_not_stressed <- go_process_summary_pdha1 %>%
  filter(stressed_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 8) # Adjust y_position as needed

single_gene_plot_pdha1 <- ggplot(go_process_summary_pdha1, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression (%)", 
       title = paste("PDHA1 Expression Compared to SHAM"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 200, by = 20), limits = c(70, 200)) +
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

single_gene_plot_pdha1




# # Create plot for DLAT
# ################################################################################################
# # Load data
# # go_process_summary_dlat <- read_excel("Data/Stress_Status/go_process_summary_dlat.xlsx")
# go_process_summary_dlat <- read_excel("Data/Stress_Status/go_process_summary_t_test_dlat.xlsx")
# 
# 
# # Setting the Timepoint as a factor to ensure order in the plot
# go_process_summary_dlat$Timepoint <- factor(go_process_summary_dlat$Timepoint,
#                                        levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
# 
# # Order the Stress_Status levels
# go_process_summary_dlat$Stress_Status <- factor(go_process_summary_dlat$Stress_Status,
#                                            levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))
# 
# 
# 
# # Modify specific timepoints if necessary
# significant_labels <- go_process_summary_dlat %>%
#   filter(significance_label != "ns") %>%
#   mutate(y_position = percentage_to_sham - 2)  # General adjustment
# 
# 
# # Filter just for significant labels; assume y_position is aligned
# stressed_not_stressed <- go_process_summary_dlat %>%
#   filter(stressed_not_stressed != "ns") %>%
#   mutate(y_position = percentage_to_sham + 2.5) # Adjust y_position as needed
# 
# 
# single_gene_plot_dlat <- ggplot(go_process_summary_dlat, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
#   geom_line(linewidth = 1.5) +
#   geom_point(size = 3) +
#   labs(x = NULL,
#        y = "Expression (%)",
#        title = paste("DLAT Expression Compared to SHAM"),
#        color = "Stress status") +  # Grouping title
#   theme_classic() +
#   scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
#   scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
#   scale_y_continuous(breaks = seq(0, 130, by = 20), limits = c(70, 130)) +
#   theme(axis.text.x = element_text(hjust = 0.5, size = 20, color = "black", margin = margin(t = 15)),
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
# single_gene_plot_dlat



# Create plot for DLD
#################################################################################################
# Load data
# go_process_summary_dld <- read_excel("Data/Stress_Status/go_process_summary_dld.xlsx")
go_process_summary_dld <- read_excel("Data/Stress_Status/go_process_summary_t_test_dld.xlsx")

# Setting the Timepoint as a factor to ensure order in the plot
go_process_summary_dld$Timepoint <- factor(go_process_summary_dld$Timepoint, 
                                           levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Order the Stress_Status levels
go_process_summary_dld$Stress_Status <- factor(go_process_summary_dld$Stress_Status,
                                               levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))

# Modify specific timepoints if necessary
significant_labels <- go_process_summary_dld %>%
  filter(significance_label != "ns") %>%
  mutate(y_position = percentage_to_sham - 0.7)  # General adjustment

# Filter just for significant labels; assume y_position is aligned
stressed_not_stressed <- go_process_summary_dld %>%
  filter(stressed_not_stressed != "ns") %>%
  mutate(y_position = percentage_to_sham + 2.5) # Adjust y_position as needed


single_gene_plot_dld <- ggplot(go_process_summary_dld, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +  
  labs(x = NULL, 
       y = "Expression (%)", 
       title = paste("DLD Expression Compared to SHAM"),
       color = "Stress status") +  # Grouping title
  theme_classic() +
  scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
  scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
  scale_y_continuous(breaks = seq(0, 160, by = 20), limits = c(60, 160)) +
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

single_gene_plot_dld



# Create plot for PDK1 
################################################################################################
# Load data
# go_process_summary_pdk1 <- read_excel("Data/Stress_Status/go_process_summary_pdk1.xlsx")
# go_process_summary_pdk1 <- read_excel("Data/Stress_Status/go_process_summary_t_test_pdk1.xlsx")
# 
# 
# 
# # Setting the Timepoint as a factor to ensure order in the plot
# go_process_summary_pdk1$Timepoint <- factor(go_process_summary_pdk1$Timepoint, 
#                                        levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
# 
# # Order the Stress_Status levels
# go_process_summary_pdk1$Stress_Status <- factor(go_process_summary_pdk1$Stress_Status,
#                                            levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))
# 
# 
# # Modify specific timepoints if necessary
# significant_labels <- go_process_summary_pdk1 %>%
#   filter(significance_label != "ns") %>%
#   mutate(y_position = ifelse(Timepoint %in% c("6 Hours") & Stress_Status %in% c("Not stressed CM"), 
#                              percentage_to_sham + 2, 
#          percentage_to_sham - 2))  # General adjustment
# 
# 
# # Filter just for significant labels; assume y_position is aligned
# stressed_not_stressed <- go_process_summary_pdk1 %>%
#   filter(stressed_not_stressed != "ns") %>%
#   mutate(y_position = percentage_to_sham + 4) # Adjust y_position as needed
# 
# 
# single_gene_plot_pdk1 <- ggplot(go_process_summary_pdk1, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
#   geom_line(linewidth = 1.5) +  
#   geom_point(size = 3) +  
#   labs(x = NULL, 
#        y = "Expression (%)", 
#        title = paste("PDK1 Expression Compared to SHAM"),
#        color = "Stress status") +  # Grouping title
#   theme_classic() +
#   scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
#   scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
#   scale_y_continuous(breaks = seq(0, 140, by = 20), limits = c(50, 140)) +
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
# single_gene_plot_pdk1



# # Create plot for PDK2 
# ################################################################################################
# # Load data
# # go_process_summary_pdk2 <- read_excel("Data/Stress_Status/go_process_summary_pdk2.xlsx")
# go_process_summary_pdk2 <- read_excel("Data/Stress_Status/go_process_summary_t_test_pdk2.xlsx")
# 
# 
# # Setting the Timepoint as a factor to ensure order in the plot
# go_process_summary_pdk2$Timepoint <- factor(go_process_summary_pdk2$Timepoint, 
#                                        levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
# 
# # Order the Stress_Status levels
# go_process_summary_pdk2$Stress_Status <- factor(go_process_summary_pdk2$Stress_Status,
#                                            levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))
# 
# 
# # Modify specific timepoints if necessary
# significant_labels <- go_process_summary_pdk2 %>%
#   filter(significance_label != "ns") %>%
#   mutate(y_position = percentage_to_sham - 3)  # General adjustment
# 
# 
# # Filter just for significant labels; assume y_position is aligned
# stressed_not_stressed <- go_process_summary_pdk2 %>%
#   filter(stressed_not_stressed != "ns") %>%
#   mutate(y_position = percentage_to_sham + 3) # Adjust y_position as needed
# 
# 
# single_gene_plot_pdk2 <- ggplot(go_process_summary_pdk2, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
#   geom_line(linewidth = 1.5) +  
#   geom_point(size = 3) +  
#   labs(x = NULL, 
#        y = "Expression (%)", 
#        title = paste("PDK2 Expression Compared to SHAM"),
#        color = "Stress status") +  # Grouping title
#   theme_classic() +
#   scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
#   scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
#   scale_y_continuous(breaks = seq(0, 140, by = 20), limits = c(40, 140)) +
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
# single_gene_plot_pdk2


# Create plot for PDHX 
################################################################################################
# Load data
# go_process_summary_pdhx <- read_excel("Data/Stress_Status/go_process_summary_pdhx.xlsx")
# go_process_summary_pdhx <- read_excel("Data/Stress_Status/go_process_summary_t_test_pdhx.xlsx")
# 
# 
# # Setting the Timepoint as a factor to ensure order in the plot
# go_process_summary_pdhx$Timepoint <- factor(go_process_summary_pdhx$Timepoint, 
#                                        levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
# 
# # Order the Stress_Status levels
# go_process_summary_pdhx$Stress_Status <- factor(go_process_summary_pdhx$Stress_Status,
#                                            levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))
# 
# 
# # Modify specific timepoints if necessary
# significant_labels <- go_process_summary_pdhx %>%
#   filter(significance_label != "ns") %>%
#   mutate(y_position = ifelse(Timepoint %in% c("1 Day"), percentage_to_sham -2, percentage_to_sham - 1.5))  # General adjustment
# 
# 
# # Filter just for significant labels; assume y_position is aligned
# stressed_not_stressed <- go_process_summary_pdhx %>%
#   filter(stressed_not_stressed != "ns") %>%
#   mutate(y_position = percentage_to_sham + 3) # Adjust y_position as needed
# 
# 
# single_gene_plot_pdhx <- ggplot(go_process_summary_pdhx, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
#   geom_line(linewidth = 1.5) +  
#   geom_point(size = 3) +  
#   labs(x = NULL, 
#        y = "Expression (%)", 
#        title = paste("PDHX Expression Compared to SHAM"),
#        color = "Stress status") +  # Grouping title
#   theme_classic() +
#   scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
#   scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
#   scale_y_continuous(breaks = seq(0, 140, by = 20), limits = c(80, 140)) +
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
# single_gene_plot_pdhx


# Combine plots:
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/pdh_multiple_panel_figure_combined.R")