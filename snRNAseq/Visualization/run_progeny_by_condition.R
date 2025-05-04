# Run PROGENy and compare ORAB vs SHAM
# The PROGENy scores are weighted scores based on predefined weights of genes. Each cell gets a score that is then summarized. 
# The scores are compared to the means of the whole dataset. 
# PROGENy scores above 0 are increased compared to the mean of the entire dataset, and below 0 are lower than the mean. 


# Source functions
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/load_libraries_progeny_analysis.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/process_data_for_progeny_by_condition.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_anova_tukey_for_progeny_scores_by_condition.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_t_test_for_progeny_scores_by_condition.R")

file.edit("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_anova_tukey_for_progeny_scores_by_condition.R")

file.edit("C:/Users/siljeew/Master_project/snRNAseq/Analysis/process_data_for_progeny_by_condition.R")

# Load libraries
load_libraries()


# Define parameters and prepare data
#########################################################################################
# Define pathways
pathways_of_interest <- c("p53", "TNFa", "NFkB", "JAK-STAT", "EGFR", "Androgen", "Hypoxia")

# Prepare paths for saving statistical test results
output_dir_plot <- "C:/Users/siljeew/Master_project/snRNAseq/Plots/PROGENy"
output_dir_data <- "C:/Users/siljeew/Master_project/snRNAseq/Data"
anova_xlsx_path <- file.path(output_dir_data, "anova_progeny_results_condition.xlsx")
tukey_xlsx_path <- file.path(output_dir_data, "tukey_progeny_results_condition.xlsx")


# Load data
mouse_vcm_all_time_points_with_stress_status <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_vcm_all_time_points_with_stress_status.Rds")

# Define the data layers in a named list
data_layers <- list(
  "6 Hours" = "data.mouse_6h_full_merge_project",
  "12 Hours" = "data.mouse_12h_full_merge_project", 
  "1 Day" = "data.mouse_1d_full_merge_project", 
  "3 Days" = "data.Full Merge Project 3d",
  "1 Week" = "data.Full Merge Project 1w",
  "3 Weeks" = "data.Full Merge Project 3w"
)

# head(mouse_vcm_all_time_points_with_stress_status)

# Process data and store in data frame
##############################################################################
conditions_df <- data.frame(Cell = colnames(mouse_vcm_all_time_points_with_stress_status),
                            orig.ident = mouse_vcm_all_time_points_with_stress_status$orig.ident,
                            Timepoint = mouse_vcm_all_time_points_with_stress_status$Timepoint,
                            Condition = mouse_vcm_all_time_points_with_stress_status$Group)


# Process each timepoint
progeny_results <- lapply(names(data_layers), function(tp) {
  process_timepoint(tp, data_layers[[tp]])}) %>%
  bind_rows() 

# Remove NaNs
progeny_results_filtered <- progeny_results %>%
  filter(Activity != "NaN")

write.xlsx(progeny_results_filtered, "C:/Users/siljeew/Master_project/snRNAseq/Data/progeny_results_by_condition.xlsx")


progeny_results_filtered$Activity <- as.numeric(progeny_results_filtered$Activity)

# Set the order as factor to define the diff correctly. This will set the diff as ORAB - SHAM. 
progeny_results_filtered$Condition <- factor(progeny_results_filtered$Condition, levels = c("SHAM", "ORAB"))

head(progeny_results_filtered)
View(progeny_results_filtered)


# Perform ANOVA+TukeyHSD test
################################################################################
# Perform the test
tukey_results_df <- perform_anova_tukey_timepoint_focus(
  data = progeny_results_filtered,
  group_col = "Condition",
  activity_col = "Activity",
  time_col = "Timepoint",
  pathways_of_interest = pathways_of_interest
)

# view(tukey_results_df)
tukey_results_df <- tukey_results_df %>%
  rename(p_adjusted = `p adj`)

# Add column for significance label
# tukey_results_df <- tukey_results_df %>%
#   mutate(significance_label = ifelse(p_adjusted < 0.05, "*",
#                                      ifelse(p_adjusted < 0.01, "**",
#                                             ifelse(p_adjusted < 0.001, "***",
#                                                    "ns"))))

view(tukey_results_df)
head(tukey_results_df)

# Write to a single Excel file with all data combined
write.xlsx(tukey_results_df, "C:/Users/siljeew/Master_project/snRNAseq/Data/tukey_results_progeny_by_condition.xlsx")



# Perform pairwise comparison
##########################################################################
# Define groups. Diff = group2 - group1. 
# group1 <- "SHAM"
# group2 <- "ORAB"
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Analysis/perform_t_test_for_progeny_scores_by_condition.R")


# Call the pairwise t-test function
t_test_results_df <- perform_pairwise_t_test(
  data = progeny_results,
  group_col = "Condition",
  activity_col = "Activity",
  time_col = "Timepoint",
  pathways_of_interest = pathways_of_interest
)

# View the results
print(t_test_results_df)


# Add significance label
t_test_results_df <- t_test_results_df %>%
  mutate(p_adjusted_label_bh = ifelse(p_adjusted_bh < 0.001, "***",
                                      ifelse(p_adjusted_bh < 0.01, "**",
                                             ifelse(p_adjusted_bh < 0.05, "*",
                                                    "ns"))))

# Save the results
write_xlsx(t_test_results_df, "C:/Users/siljeew/Master_project/snRNAseq/Data/t_test_result_progeny_by_condition.xlsx")


# Create plot
###########################################################################################################################
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/create_plot_for_combined_progeny_scores_by_condition.R")



