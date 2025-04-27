# PROGENy

library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tidyverse)
library(tibble)
library(writexl)

# Inspect layers 
Layers(mouse_vcm_all_time_points_with_stress_status)


# Extract Data
data_matrix_3w <- GetAssayData(object = mouse_vcm_all_time_points_with_stress_status, layer = "data.Full Merge Project 3w", assay = "RNA")
dense_matrix <- as.matrix(data_matrix_3w)

# Run Progeny Analysis
pathway_activity_scores <- progeny(dense_matrix, scale = FALSE, organism = "Mouse", top = 500)
scaled_pathway_scores <- scale(pathway_activity_scores)
progeny_scores_df <- as.data.frame(t(scaled_pathway_scores)) %>%
  tibble::rownames_to_column("Cell") %>%
  tidyr::gather(Pathway, Activity, -Cell)
progeny_scores_df <- progeny_scores_df %>%
  rename(Pathway = Cell, Cell = Pathway)

# Print Data Frame
View(progeny_scores_df)

# Extract meta data information
head(mouse_vcm_all_time_points_with_stress_status@meta.data$Stress_Status)
conditions_df <- data.frame(Cell = colnames(mouse_vcm_all_time_points_with_stress_status),
                            Stress_Status = mouse_vcm_all_time_points_with_stress_status$Stress_Status,
                            Timepoint = mouse_vcm_all_time_points_with_stress_status$Timepoint)

# Ensure pathway activity data frame includes a "Cell" column
progeny_scores_df <- progeny_scores_df %>%
  inner_join(conditions_df, by = "Cell")  # Join to add condition info to each cell

#View(progeny_scores_df)




# Group and summarise progeny scores
######################################################################
# Group by Stress_Status instead of cell type
summarized_progeny_scores <- progeny_scores_df %>%
  group_by(Pathway, Stress_Status, Timepoint) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

# Prepare the data frame for heatmap plotting
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%
  tidyr::spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
summarized_progeny_scores_df[is.na(summarized_progeny_scores_df)] <- 0

# Verify the modified data frame
print(summarized_progeny_scores_df)
view(summarized_progeny_scores_df)

progeny_scores_3w <- summarized_progeny_scores_df
# Convert row names into a column named "Identifier"
progeny_scores_3w <- summarized_progeny_scores_df %>%
  tibble::rownames_to_column(var = "Stress_Status")

# Verify the changes
head(progeny_scores_3w)










# Compare stressed and SHAM
##############################################################################
# Subset data
head(progeny_scores_df)
sham_means <- progeny_scores_df %>%
  filter(Stress_Status == "SHAM - CM") %>%
  group_by(Pathway, Timepoint) %>%
  summarise(mean_activity = mean(Activity, na.rm = TRUE))
view(sham_means)

# not_stressed_means <- progeny_scores_df %>%
#   filter(Stress_Status == "Not stressed CM") %>%
#   group_by(Pathway, Timepoint) %>%
#   summarise(mean_activity = mean(Activity, na.rm = TRUE))
# view(not_stressed_means)

stressed_means <- progeny_scores_df %>%
  filter(Stress_Status == "Stressed CM") %>%
  group_by(Pathway, Timepoint) %>%
  summarise(mean_activity = mean(Activity, na.rm = TRUE))
view(stressed_means)



# Join the datasets on Pathway and Timepoint
activity_comparison <- left_join(sham_means, stressed_means, 
                                 by = c("Pathway", "Timepoint"),
                                 suffix = c("_sham", "_stressed"))
head(activity_comparison)

# Calculate Differences
activity_comparison <- activity_comparison %>%
  mutate(difference = mean_activity_stressed - mean_activity_sham)

# Print differences
View(activity_comparison)

# Remove NA values and zero differences
activity_filtered <- activity_comparison %>%
  filter(!is.na(difference), difference != 0)

# # Determine the maximum absolute difference
# max_diff <- max(abs(activity_filtered$difference), na.rm = TRUE)

# Create the plot
activity_filtered_3w <- activity_filtered

# Import P-adjusted label
activity_filtered_3w <- activity_filtered_3w %>%
  left_join(updated_results_df %>% select(Pathway, Timepoint, p_adjusted_label), 
            by = c("Pathway", "Timepoint"))

# Remove NAs
p_adj_label_3w <- activity_filtered_3w %>%
  filter(p_adjusted_label != "ns") %>%
  mutate(y_position = ifelse(Pathway == c("JAK-STAT"), 
                             difference -0.02, 
                             difference + 0.01))

# Create the plot
progeny_3w <- ggplot(activity_filtered_3w, aes(x = Pathway, y = difference, fill = difference > 0)) +
  geom_bar(stat = "identity") +
  labs(x = "", 
       y = "Difference in activity score",
       title = "SHAM vs Stressed CM at 3 Weeks",
       fill = "Condition") +
  scale_fill_manual(
    values = c("TRUE" = "coral", "FALSE" = "grey22"),
    labels = c("Stressed CM", "SHAM CM"),
    breaks = c(TRUE, FALSE)
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +  
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0.5, size = 16, color = "black", margin = margin(t = 15)), 
        axis.text.y = element_text(size = 26, colour = "black"), 
        axis.title.y = element_text(size = 26),
        axis.title.x = element_blank(), 
        plot.title = element_text(size = 26, margin = margin(b = 20)),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 26),
        panel.border = element_blank(),
        axis.line.y.right = element_blank()) +
  geom_text(data = p_adj_label_3w, aes(x = Pathway, y = y_position, label = p_adjusted_label), 
            size = 10, color = "black") +
  ylim(-0.2, 0.25)

progeny_3w


# Perform t-test
################################################################################
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_t_test_for_progeny_scores.R")
progeny_scores_t_test <- read.xlsx("C:/Users/siljeew/Master_project/snRNAseq/Data/progeny_scores_t_test.xlsx")
View(progeny_scores_t_test)

# Subset data
head(progeny_scores_df)

# # Define pathway of interest
# pathway_of_interest <- "p53"
# pathway_of_interest <- "TNFa"
# pathway_of_interest <- "NFkB"
# pathway_of_interest <- "JAK-STAT"
# pathway_of_interest <- "EGFR"
# pathway_of_interest <- "Androgen"
# 
# 
# # Define the groups to compare
# sham_group <- "SHAM - CM"
# stressed_group <- "Stressed CM"
# 
# 
# # Filter data for the specific pathway and the relevant groups
# data_for_t_test <- progeny_scores_df %>%
#   filter(Stress_Status %in% c(sham_group, stressed_group)) %>%  # Filter the Stress_Status
#   filter(Pathway == pathway_of_interest) %>%  # Filter the pathway of interest
#   group_by(Pathway) %>%
#   filter(n_distinct(Stress_Status) == 2) %>%
#   ungroup()
# 
# 
# head(data_for_t_test)
# 
# # Perform t-test for the specific Pathway
# test_result <- data_for_t_test %>%
#   summarise(t_test = list(t.test(Activity ~ Stress_Status, data = .))) %>%
#   ungroup()
# 
# head(test_result)
# #View(test_result)
# 
# # Extract and view the p-value for the specific t-test
# test_pvalue <- sapply(test_result$t_test, function(x) x$p.value)
# print(test_pvalue)








