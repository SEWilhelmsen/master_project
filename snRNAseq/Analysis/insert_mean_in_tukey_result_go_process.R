# Combine Tukey results with mean expression
###############################################################################
anova_tukey_process_result <- anova_tukey_process_result %>%
  dplyr::select(-Timepoint)

anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group1_mean = NA,
         group2_mean = NA)

  

# ORAB 6 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(process == go_process_of_interest, time_point == "6 Hours", condition == "ORAB") %>%
  pull(mean_expression) 


# Update "group1_mean"
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group1_mean = ifelse(process == go_process_of_interest & group1_timepoint == "6 Hours" & group1 == "ORAB", value, group1_mean))

# Update "group2_mean"
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group2_mean = ifelse(process == go_process_of_interest & group2_timepoint == "6 Hours" & group2 == "ORAB", value, group2_mean))


# ORAB 12 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(process == go_process_of_interest, time_point == "12 Hours", condition == "ORAB") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean"
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group1_mean = ifelse(process == go_process_of_interest & group1_timepoint == "12 Hours" & group1 == "ORAB", value, group1_mean))

# Update "group1_mean"
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group2_mean = ifelse(process == go_process_of_interest & group2_timepoint == "12 Hours" & group2 == "ORAB", value, group2_mean))


# ORAB 1 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(process == go_process_of_interest, time_point == "1 Day", condition == "ORAB") %>%
  pull(mean_expression) 

# Update "group1_mean"
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group1_mean = ifelse(process == go_process_of_interest & group1_timepoint == "1 Day" & group1 == "ORAB", value, group1_mean))


# Update "group1_mean"
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group2_mean = ifelse(process == go_process_of_interest & group2_timepoint == "1 Day" & group2 == "ORAB", value, group2_mean))



# ORAB 3 Days
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(process == go_process_of_interest, time_point == "3 Days", condition == "ORAB") %>%
  pull(mean_expression) 

# Update "group1_mean"
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group1_mean = ifelse(process == go_process_of_interest & group1_timepoint == "3 Days" & group1 == "ORAB", value, group1_mean))

# Update "group2_mean"
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group2_mean = ifelse(process == go_process_of_interest & group2_timepoint == "3 Days" & group2 == "ORAB", value, group2_mean))



# ORAB 1 Week
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(process == go_process_of_interest, time_point == "1 Week", condition == "ORAB") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean"
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group1_mean = ifelse(process == go_process_of_interest & group1_timepoint == "1 Week" & group1 == "ORAB", value, group1_mean))

# Update "group2_mean"
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group2_mean = ifelse(process == go_process_of_interest & group2_timepoint == "1 Week" & group2 == "ORAB", value, group2_mean))



# Stressed 3 Weeks
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(process == go_process_of_interest, time_point == "3 Weeks", condition == "ORAB") %>%
  pull(mean_expression) 

# Update "group1_mean" 
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group1_mean = ifelse(process == go_process_of_interest & group1_timepoint == "3 Weeks" & group1 == "ORAB", value, group1_mean))

# Update "group2_mean" 
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group2_mean = ifelse(process == go_process_of_interest & group2_timepoint == "3 Weeks" & group2 == "ORAB", value, group2_mean))


view(anova_tukey_process_result)




# SHAM
################################################################################
# 6 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(process == go_process_of_interest, time_point == "6 Hours", condition == "SHAM") %>%
  pull(mean_expression) 

# Update "group1_mean" 
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group1_mean = ifelse(process == go_process_of_interest & group1_timepoint == "6 Hours" & group1 == "SHAM", value, group1_mean))

# Update "group2_mean" 
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group2_mean = ifelse(process == go_process_of_interest & group2_timepoint == "6 Hours" & group2 == "SHAM", value, group2_mean))



# SHAM 12 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(process == go_process_of_interest, time_point == "12 Hours", condition == "SHAM") %>%
  pull(mean_expression) 

# Update "group1_mean"
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group1_mean = ifelse(process == go_process_of_interest & group1_timepoint == "12 Hours" & group1 == "SHAM", value, group1_mean))

# Update "group2_mean"
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group2_mean = ifelse(process == go_process_of_interest & group2_timepoint == "12 Hours" & group2 == "SHAM", value, group2_mean))



# SHAM 1 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(process == go_process_of_interest, time_point == "1 Day", condition == "SHAM") %>%
  pull(mean_expression) 


# Update "group1_mean" 
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group1_mean = ifelse(process == go_process_of_interest & group1_timepoint == "1 Day" & group1 == "SHAM", value, group1_mean))


# Update "group2_mean" 
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group2_mean = ifelse(process == go_process_of_interest & group2_timepoint == "1 Day" & group2 == "SHAM", value, group2_mean))



# SHAM 3 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(process == go_process_of_interest, time_point == "3 Days", condition == "SHAM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean" 
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group1_mean = ifelse(process == go_process_of_interest & group1_timepoint == "3 Days" & group1 == "SHAM", value, group1_mean))

# Update "group2_mean" 
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group2_mean = ifelse(process == go_process_of_interest & group2_timepoint == "3 Days" & group2 == "SHAM", value, group2_mean))


# SHAM 1 Week
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(process == go_process_of_interest, time_point == "1 Week", condition == "SHAM") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean" 
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group1_mean = ifelse(process == go_process_of_interest & group1_timepoint == "1 Week" & group1 == "SHAM", value, group1_mean))

# Update "group2_mean" 
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group2_mean = ifelse(process == go_process_of_interest & group2_timepoint == "1 Week" & group2 == "SHAM", value, group2_mean))


# SHAM 3 Weeks
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(process == go_process_of_interest, time_point == "3 Weeks", condition == "SHAM") %>%
  pull(mean_expression) 

# Update "group1_mean" 
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group1_mean = ifelse(process == go_process_of_interest & group1_timepoint == "3 Weeks" & group1 == "SHAM", value, group1_mean))

# Update "group2_mean" 
anova_tukey_process_result <- anova_tukey_process_result %>%
  mutate(group2_mean = ifelse(process == go_process_of_interest & group2_timepoint == "3 Weeks" & group2 == "SHAM", value, group2_mean))


view(anova_tukey_process_result)


