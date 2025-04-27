# PROGENy

library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tidyverse)
library(writexl)
library(tibble)
library(openxlsx)


# Load data
mouse_vcm_all_time_points_with_stress_status <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_vcm_all_time_points_with_stress_status.Rds")

# Inspect layers 
Layers(mouse_vcm_all_time_points_with_stress_status)

# Extract Data
data_matrix_6h <- GetAssayData(object = mouse_vcm_all_time_points_with_stress_status, layer = "data.mouse_6h_full_merge_project", assay = "RNA")
dense_matrix <- as.matrix(data_matrix_6h)

# Run Progeny Analysis
pathway_activity_scores <- progeny(dense_matrix, scale = FALSE, organism = "Mouse", top = 500)
scaled_pathway_scores <- scale(pathway_activity_scores)
progeny_scores_df <- as.data.frame(t(scaled_pathway_scores)) %>%
  tibble::rownames_to_column("Cell") %>%
  tidyr::gather(Pathway, Activity, -Cell)
progeny_scores_df <- progeny_scores_df %>%
  rename(Pathway = Cell, Cell = Pathway)
  
View(progeny_scores_df)

# Print Data Frame
print(progeny_scores_df)

# Extract condition information from Seurat metadata
head(mouse_vcm_all_time_points_with_stress_status@meta.data$Stress_Status)
conditions_df <- data.frame(Cell = colnames(mouse_vcm_all_time_points_with_stress_status),
                            Stress_Status = mouse_vcm_all_time_points_with_stress_status$Stress_Status,
                            Timepoint = mouse_vcm_all_time_points_with_stress_status$Timepoint)

# Ensure pathway activity data frame includes a "Cell" column
progeny_scores_df <- progeny_scores_df %>%
  inner_join(conditions_df, by = "Cell")  # Join to add condition info to each cell

head(conditions_df)
View(progeny_scores_df)


# Prepare data for heatmap
#####################################################################
# Group by Stress_Status instead of cell type
summarized_progeny_scores <- progeny_scores_df %>%
  group_by(Pathway, Stress_Status, Timepoint) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

View(summarized_progeny_scores)


# Prepare the data frame for heatmap plotting
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%
  tidyr::spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
summarized_progeny_scores_df[is.na(summarized_progeny_scores_df)] <- 0

# Verify the modified data frame
print(summarized_progeny_scores_df)
view(summarized_progeny_scores_df)

progeny_scores_6h <- summarized_progeny_scores_df

# Convert row names into a column named "Identifier"
progeny_scores_6h <- summarized_progeny_scores_df %>%
  tibble::rownames_to_column(var = "Stress_Status")

# Verify the changes
head(progeny_scores_6h)
view(progeny_scores_6h)



# Store results in a data frame 
############################################################################
# Prepare the data frame for heatmap plotting with Timepoint
summarized_progeny_scores_for_xlsx <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%
  tidyr::spread(Pathway, avg) %>%
  # Ensure Stress_Status and Timepoint are still columns in the final data frame
  dplyr::ungroup() %>%  # Ungroup if further processing without grouping is required
  data.frame(stringsAsFactors = FALSE)  

# Verify the modified data frame
print(summarized_progeny_scores_for_xlsx)
View(summarized_progeny_scores_for_xlsx)
# Save to xlsx
write_xlsx(summarized_progeny_scores_for_xlsx, path = "C:/Users/siljeew/Master_project/snRNAseq/Data/progeny_scores.xlsx")




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
# View(activity_comparison)

# Remove NA values and zero differences
activity_filtered <- activity_comparison %>%
  filter(!is.na(difference), difference != 0)

# Determine the maximum absolute difference
max_diff <- max(abs(activity_filtered$difference), na.rm = TRUE)
head(activity_filtered)

activity_filtered_6h <- activity_filtered

# Statistical test
################################################################################
# Perform t-test for means
# source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_t_test_for_progeny_scores.R")
# progeny_scores_t_test <- read.xlsx("C:/Users/siljeew/Master_project/snRNAseq/Data/progeny_scores_t_test.xlsx")
# View(progeny_scores_t_test)


# Import P-adjusted label
activity_filtered_6h <- activity_filtered_6h %>%
  left_join(updated_results_df %>% select(Pathway, Timepoint, p_adjusted_label), 
            by = c("Pathway", "Timepoint"))

# Remove NAs
p_adj_label_6h <- activity_filtered_6h %>%
  filter(p_adjusted_label != "ns") %>%
  mutate(y_position = difference + 0.5)
  
print(p_adj_label_6h)  

# Create plot
progeny_6h <- ggplot(activity_filtered_6h, aes(x = Pathway, y = difference, fill = difference > 0)) +
  geom_bar(stat = "identity") +
  labs(x = "", 
       y = "Difference in activity score",
       title = "SHAM vs Stressed CM at 6 Hours") +
  scale_fill_manual(values = c("TRUE" = "coral", "FALSE" = "grey22")) +
  guides(fill = "none") + 
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
  geom_text(data = p_adj_label_6h, aes(x = Pathway, y = y_position, label = p_adjusted_label), 
            size = 10, vjust = -0.5, color = "black") +
  ylim(-0.2, 0.25)

progeny_6h






# With only one pathway, but all stress groups
###############################################################################
# Specify the pathway you want to visualize
single_pathway <- "EGFR"  # Replace with the actual pathway name

# Filter out NA values from the relevant column(s)
filtered_data_no_na <- progeny_scores_df %>%
  filter(!is.na(Activity))

# Filter the data to include only the specific pathway
filtered_single_pathway <- filtered_data_no_na %>%
  filter(Pathway == single_pathway)

View(filtered_data_no_na)

# Create the violin plot for the single pathway
ggplot(filtered_single_pathway, aes(x = Pathway, y = Activity, fill = Stress_Status)) +
  geom_violin(trim = TRUE) +
  scale_fill_manual(values = c("SHAM - CM" = "grey22", "Stressed CM" = "coral", "Not stressed CM" = "red")) +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_blank()
  ) +
  labs(y = "Activity Score") +
  facet_wrap(~ Timepoint, scales = "free_x")


