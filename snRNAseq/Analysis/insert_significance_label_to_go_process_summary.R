# Insert significance labels and adjusted p-value to go_process_summary from anova_tukey_process_results
# Create new columns in go_process_summary, pull values and insert in new columns.
# Silje Wilhelmsen

# Create column for significance labels and set as factors
############################################################################
go_process_summary <- go_process_summary %>%
  mutate(significance_label = NA, 
         p_adj = NA)

View(go_process_summary)

go_process_summary$time_point <- factor(go_process_summary$time_point)
go_process_summary$p_adj <- factor(go_process_summary$p_adj)



# 6 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_process_result %>%
  filter(process == go_process_of_interest, group1_timepoint == "6 Hours", group2_timepoint == "6 Hours", group1 != group2) %>%
  pull(significance_label) 

value <- anova_tukey_process_result %>%
  filter(process == go_process_of_interest, group1_timepoint == "6 Hours", group2_timepoint == "6 Hours", group1 != group2) %>%
  pull("p adj") 

# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(process == go_process_of_interest & time_point == "6 Hours", padj_label, significance_label)) %>%
  mutate(p_adj = ifelse(process == go_process_of_interest & time_point == "6 Hours", value, p_adj))



# 12 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_process_result %>%
  filter(process == go_process_of_interest, group1_timepoint == "12 Hours", group2_timepoint == "12 Hours", group1 != group2) %>%
  pull(significance_label) 

value <- anova_tukey_process_result %>%
  filter(process == go_process_of_interest, group1_timepoint == "12 Hours", group2_timepoint == "12 Hours", group1 != group2) %>%
  pull("p adj") 

# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(process == go_process_of_interest & time_point == "12 Hours", padj_label, significance_label)) %>%
  mutate(p_adj = ifelse(process == go_process_of_interest & time_point == "12 Hours", value, p_adj))




# 1 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_process_result %>%
  filter(process == go_process_of_interest, group1_timepoint == "1 Day", group2_timepoint == "1 Day", group1 != group2) %>%
  pull(significance_label) 

value <- anova_tukey_process_result %>%
  filter(process == go_process_of_interest, group1_timepoint == "1 Day", group2_timepoint == "1 Day", group1 != group2) %>%
  pull("p adj") 

# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(process == go_process_of_interest & time_point == "1 Day", padj_label, significance_label)) %>%
  mutate(p_adj = ifelse(process == go_process_of_interest & time_point == "1 Day", value, p_adj))




# 3 Days
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_process_result %>%
  filter(process == go_process_of_interest, group1_timepoint == "3 Days", group2_timepoint == "3 Days", group1 != group2) %>%
  pull(significance_label) 

value <- anova_tukey_process_result %>%
  filter(process == go_process_of_interest, group1_timepoint == "3 Days", group2_timepoint == "3 Days", group1 != group2) %>%
  pull("p adj") 

# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(process == go_process_of_interest & time_point == "3 Days", padj_label, significance_label)) %>%
  mutate(p_adj = ifelse(process == go_process_of_interest & time_point == "3 Days", value, p_adj))




# 1 Week
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_process_result %>%
  filter(process == go_process_of_interest, group1_timepoint == "1 Week", group2_timepoint == "1 Week", group1 != group2) %>%
  pull(significance_label) 

value <- anova_tukey_process_result %>%
  filter(process == go_process_of_interest, group1_timepoint == "1 Week", group2_timepoint == "1 Week", group1 != group2) %>%
  pull("p adj") 

# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(process == go_process_of_interest & time_point == "1 Week", padj_label, significance_label)) %>%
  mutate(p_adj = ifelse(process == go_process_of_interest & time_point == "1 Week", value, p_adj))




# 3 Weeks
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_process_result %>%
  filter(process == go_process_of_interest, group1_timepoint == "3 Weeks", group2_timepoint == "3 Weeks", group1 != group2) %>%
  pull(significance_label) 

value <- anova_tukey_process_result %>%
  filter(process == go_process_of_interest, group1_timepoint == "3 Weeks", group2_timepoint == "3 Weeks", group1 != group2) %>%
  pull("p adj") 

# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(process == go_process_of_interest & time_point == "3 Weeks", padj_label, significance_label)) %>%
  mutate(p_adj = ifelse(process == go_process_of_interest & time_point == "3 Weeks", value, p_adj))