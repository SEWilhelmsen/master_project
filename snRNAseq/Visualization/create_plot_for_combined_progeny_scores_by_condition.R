# Create plot for PROGENy activity scores by condition
# Silje Wilhelmsen

library(ggplot2)
library(dplyr)
library(tidyr)

# Data prepared from herer:
#file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/run_progeny_by_condition.R")

# Colors
# SHAM: "grey22"
# ORAB: "coral"

# Verify that data looks correct
# View(t_test_results_df)
# head(t_test_results_df)
# head(tukey_results_df)


# Prepare data for plot
#############################################################
# Insert Benjamini-Hochberg adjusted p-value from t_test_results_df to progeny_for_plot
progeny_for_plot <- tukey_results_df %>%
  left_join(t_test_results_df %>%
              dplyr::select(Pathway, Timepoint, p_adjusted_bh, p_adjusted_label_bh),
            by = c("Pathway", "Timepoint"))

View(progeny_for_plot)
head(progeny_for_plot)

# # Create object of labels
# p_adj_label <- progeny_for_plot %>%
#   filter(p_adjusted_label_bh != "ns") %>%
#   mutate(x_position = ifelse(diff < 0 & Timepoint %in% c("12 Hours", "3 Days"), diff - 0.05,
#                              ifelse(diff < 0 & Timepoint %in% c("3 Weeks"), diff - 0.03,
#                                     ifelse(diff > 0 & Timepoint %in% c("1 Day"), diff +0.08,
#                                            diff + 0.02))))

# Or all stars on the same side
p_adj_label <- progeny_for_plot %>%
  filter(p_adjusted_label_bh != "ns") %>%
  mutate(x_position = ifelse(diff > 0 & Timepoint %in% c("12 Hours", "3 Days"), upr + 0.05,
                             ifelse(diff > 0 & Timepoint %in% c("3 Weeks"), upr + 0.03,
                                    ifelse(diff < 0 & Timepoint %in% c("1 Day"), lwr - 0.06,
                                           lwr - 0.03))))

# Set timepoint as factor for correct order
progeny_for_plot$Timepoint <- factor(progeny_for_plot$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
p_adj_label$Timepoint <- factor(p_adj_label$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))



# Create plot
##############################################################
# Create plot
progeny_plot <- ggplot(progeny_for_plot, aes(x = diff, y = Pathway, fill = diff > 0)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(y = Pathway, xmin = lwr, xmax = upr),
                width = 0.5, linewidth = 1, 
                position = position_dodge(0.9)) +  # Adjusted to match the bars
  labs(x = "Difference in activity score",
       y = "",
       title = "SHAM vs ORAB") +
  facet_wrap(~ Timepoint, ncol = 2) +
  scale_fill_manual(values = c("TRUE" = "coral", "FALSE" = "grey22")) +
  guides(fill = "none") +
  theme_classic() +
  theme(plot.title = element_text(size = 30, margin = margin(b = 15, t = 10), hjust = 0.5),
        strip.text = element_text(size = 30, color = "black"),  # Adjust timepoint title
        panel.border = element_blank(),
        axis.title.x = element_text(size = 26, color = "black", margin = margin(t = 15)),
        axis.text.x = element_text(hjust = 0.5, size = 20, color = "black", margin = margin(t = 5)),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 28, colour = "black", margin = margin(r = 5)),
        axis.line.y.right = element_blank(),
        strip.background = element_blank()) +  # Remove strip background
  geom_text(data = p_adj_label, aes(y = Pathway, x = x_position, label = p_adjusted_label_bh), 
            size = 18, vjust = 0.8, color = "black") +
  xlim(-0.35, 0.3)

progeny_plot

ggsave(file.path(output_dir_plot, "combined_progeny_condition.png"), plot = progeny_plot, width = 15, height = 18)
ggsave(file.path(output_dir_plot, "combined_progeny_condition.pdf"), plot = progeny_plot, width = 15, height = 18)



# Create a legend using a separate data frame
legend_df <- data.frame(
  Condition = c("ORAB", "SHAM"),
  color = c("coral", "grey22"),
  diff = c(1, 1)  # Using dummy values
)

legend_df$Condition <- factor(legend_df$Condition, levels = c("SHAM", "ORAB"))

# Create a separate plot for legend
legend_plot <- ggplot(legend_df, aes(x = factor(1), y = diff, fill = Condition)) +
  geom_bar(stat = "identity") +
  labs(fill = "Condition") +
  scale_fill_manual(values = c("SHAM" = "grey22", "ORAB" = "coral"), 
                    labels = c("Higher activity score in SHAM", "Higher activity score in ORAB")) +
  theme_void() +
  theme(legend.title = element_text(size = 30),  
        legend.text = element_text(size = 30))

# Extract the legends
legend <- get_legend(legend_plot)

# Combine the main plot with the manual legend
combined_plot <- plot_grid(progeny_plot, legend, ncol = 1, rel_heights = c(8, 1))

# Display the final plot
combined_plot

ggsave(file.path(output_dir_plot, "combined_progeny_plot_legend_orab_vs_sham.png"), plot = combined_plot, width = 18, height = 20)
ggsave(file.path(output_dir_plot, "combined_progeny_plot_legend_orab_vs_sham.pdf"), plot = combined_plot, width = 18, height = 20)


head(t_test_results_df)

