# library(dplyr)
# library(ggplot2)
# library(ggpubr)
# library(broom)
# 
# #Use the results of the anova_tukey_gene_results to manually add p.adj and label to go_process_summary
# 
# 
# gene_expression_summary <- read_excel("Data/Stress_Status/go_process_summary.xlsx")
# View(gene_expression_summary)
# 
# gene_expression_summary$Stress_Status <- factor(gene_expression_summary$Stress_Status, levels = c("SHAM CM", "Not stressed CM", "Stressed CM"))
# 
# # Add a y_position column to gene_expression_summary for labels to appear above bars
# gene_expression_summary <- gene_expression_summary %>%
#   mutate(y_position = percentage_to_sham + 1) 
# View(gene_expression_summary)
# 
# # Filter labels_df to exclude 'ns' labels
# filtered_labels_df <- gene_expression_summary %>%
#   filter(label != "ns")
# View(filtered_labels_df)
# 
# 
# # Plot with added significance annotations excluding 'ns'
# single_gene_plot <- ggplot(gene_expression_summary, aes(x = Timepoint, y = percentage_to_sham, fill = Stress_Status)) +
#   geom_bar(stat = "identity", width = 0.9, position = position_dodge()) +
#   labs(x = NULL, y = "Expression compared to SHAM CM (%)", 
#        title = paste("Expression of", gene_of_interest)) +
#   theme_classic() +
#   scale_fill_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
#   scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
#   scale_y_break(c(10, 60)) +
#   scale_y_continuous(breaks = seq(0, 260, by = 20)) + 
#   theme(axis.text.x = element_text(hjust = 0.5, size = 16, color = "black"),
#         axis.text.y = element_text(size = 12, colour = "black"),
#         axis.line.y.right = element_blank()) +
#   geom_text(data = filtered_labels_df, aes(x = Timepoint, y = y_position, label = label), 
#             position = position_dodge(width = 0.6), 
#             vjust = 0, hjust = -1, size = 10)
# 
# single_gene_plot
# 
# single_gene_plot_bdh1 <- ggplot(gene_expression_summary, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
#   geom_line(size = 1.5) +  
#   geom_point(size = 3) +  
#   labs(x = NULL, y = "Expression compared to SHAM CM (%)", 
#        title = paste("Expression of", gene_of_interest, "compared to SHAM")) +
#   theme_classic() +
#   scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
#   scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
#   scale_y_continuous(breaks = seq(0, 260, by = 20)) + 
#   theme(axis.text.x = element_text(hjust = 0.5, size = 16, color = "black"),
#         axis.text.y = element_text(size = 12, colour = "black"),
#         axis.line.y.right = element_blank()) +
#   geom_text(data = filtered_labels_df, aes(x = Timepoint, y = y_position, label = label), size = 10)
# 
# single_gene_plot_bdh1
# 
# 
# 
# 
# # Include significance triangle if p.adj < 0.05 between stressed and not stressed
# #####################################################################################
# # Filter and select columns based on the specified condition
# filtered_anova_tukey_results <- anova_tukey_gene_results %>%
#   filter((group1 %in% c("Stressed CM", "Not stressed CM")) &
#            (group2 %in% c("Stressed CM", "Not stressed CM"))) %>%
#   dplyr::select(Timepoint, gene, `p adj`) %>%
#   mutate(label = case_when(
#     `p adj` < 0.05 ~ "\u25B2",  # Unicode for triangle
#     TRUE ~ "ns"
#   )) %>%
#   mutate(Stress_Status = "Stressed CM")
# 
# # Print or view the filtered data
# print(filtered_anova_tukey_results)
# 
# correct_levels <- c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")
# 
# # Filter for only significant results
# significant_results <- filtered_anova_tukey_results %>%
#   filter(label != "ns")
# print(significant_results)
# 
# # WTF er trekanten 
# single_gene_plot_bdh1 <- ggplot(gene_expression_summary, aes(x = Timepoint, y = percentage_to_sham, color = Stress_Status, group = Stress_Status)) +
#   geom_line(size = 1.5) +  
#   geom_point(size = 3) +  
#   labs(x = NULL, y = "Expression compared to SHAM CM (%)", 
#        title = paste("Expression of", gene_of_interest, "compared to SHAM")) +
#   theme_classic() +
#   scale_color_manual(values = c("SHAM CM" = "grey22", "Not stressed CM" = "darkorange", "Stressed CM" = "coral3")) +
#   scale_x_discrete(labels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")) +
#   scale_y_continuous(breaks = seq(0, 260, by = 20)) + 
#   theme(axis.text.x = element_text(hjust = 0.5, size = 16, color = "black"),
#         axis.text.y = element_text(size = 12, colour = "black"),
#         axis.line.y.right = element_blank()) +
#   geom_text(data = filtered_labels_df, aes(x = Timepoint, y = y_position, label = label), size = 10) +
#   geom_point(data = significant_results, aes(x = Timepoint), 
#              y = 260, shape = 2, size = 6, color = "coral3")
# 
# single_gene_plot_bdh1
# 
