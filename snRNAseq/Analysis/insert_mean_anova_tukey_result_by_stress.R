# Combine ANOVA and Tukey results with mean expression
###############################################################################

anova_tukey_mean <- anova_tukey_gene_results %>% 
  mutate(group1_mean_expression = NA)

# Stressed 6 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "6 Hours", Stress_Status == "Stressed CM") %>%
  pull(mean_expression) 

# Print the retrieved value
print(value)

colnames(anova_tukey_mean)

# Update "group1_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "6 Hours" & group1 == "Stressed CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "6 Hours", group1 == "Stressed CM"))



# Stressed 12 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "12 Hours", Stress_Status == "Stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "12 Hours" & group1 == "Stressed CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "12 Hours", group1 == "Stressed CM"))



# Stressed 1 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "1 Day", Stress_Status == "Stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "1 Day" & group1 == "Stressed CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "1 Day", group1 == "Stressed CM"))



# Stressed 3 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "3 Days", Stress_Status == "Stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "3 Days" & group1 == "Stressed CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "3 Days", group1 == "Stressed CM"))



# Stressed 1 Week
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "1 Week", Stress_Status == "Stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "1 Week" & group1 == "Stressed CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "1 Week", group1 == "Stressed CM"))



# Stressed 3 Weeks
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "3 Weeks", Stress_Status == "Stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "3 Weeks" & group1 == "Stressed CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "3 Weeks", group1 == "Stressed CM"))


view(anova_tukey_mean)




# Not stressed
################################################################################
# Not stressed 6 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "6 Hours", Stress_Status == "Not stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "6 Hours" & group1 == "Not stressed CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "6 Hours", group1 == "Not stressed CM"))



# Not stressed 12 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "12 Hours", Stress_Status == "Not stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "12 Hours" & group1 == "Not stressed CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "12 Hours", group1 == "Not stressed CM"))



# Not stressed 1 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "1 Day", Stress_Status == "Not stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "1 Day" & group1 == "Not stressed CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "1 Day", group1 == "Not stressed CM"))



# Not stressed 3 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "3 Days", Stress_Status == "Not stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "3 Days" & group1 == "Not stressed CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "3 Days", group1 == "Not stressed CM"))


# Not stressed 1 Week
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "1 Week", Stress_Status == "Not stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "1 Week" & group1 == "Not stressed CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "1 Week", group1 == "Not stressed CM"))


# Not Stressed 3 Weeks
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "3 Weeks", Stress_Status == "Not stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "3 Weeks" & group1 == "Not stressed CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "3 Weeks", group1 == "Not stressed CM"))


view(anova_tukey_mean)









#############################################################################
# SHAM
#############################################################################
# sham 6 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "6 Hours", Stress_Status == "SHAM CM") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "6 Hours" & group1 == "SHAM CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "6 Hours", group1 == "SHAM CM"))



# sham 12 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "12 Hours", Stress_Status == "SHAM CM") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "12 Hours" & group1 == "SHAM CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "12 Hours", group1 == "SHAM CM"))



# sham 1 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "1 Day", Stress_Status == "SHAM CM") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "1 Day" & group1 == "SHAM CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "1 Day", group1 == "SHAM CM"))



# sham 3 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "3 Days", Stress_Status == "SHAM CM") %>%
  pull(mean_expression) 

# Print the retrieved value
print(value)

colnames(anova_tukey_mean)
# Update "group1_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "3 Days" & group1 == "SHAM CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "3 Days", group1 == "SHAM CM"))



# sham 1 Week
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "1 Week", Stress_Status == "SHAM CM") %>%
  pull(mean_expression) 

# Print the retrieved value
print(value)

# Update "group1_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "1 Week" & group1 == "SHAM CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "1 Week", group1 == "SHAM CM"))



# SHAM 3 Weeks
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "3 Weeks", Stress_Status == "SHAM CM") %>%
  pull(mean_expression) 
print(value)

# Update "group1_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group1_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "3 Weeks" & group1 == "SHAM CM", value, group1_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "3 Weeks", group1 == "SHAM CM"))








#############################################################################
# Group2
#############################################################################
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = NA)
View(anova_tukey_mean)

# anova_tukey_mean <- anova_tukey_mean %>%
#   dplyr::select(-group2_mean_expression)
# View(anova_tukey_mean)


# Stressed 6 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "6 Hours", Stress_Status == "Stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "6 Hours" & group2 == "Stressed CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "6 Hours", group2 == "Stressed CM"))



# Stressed 12 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "12 Hours", Stress_Status == "Stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "12 Hours" & group2 == "Stressed CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "12 Hours", group2 == "Stressed CM"))




# Stressed 1 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "1 Day", Stress_Status == "Stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "1 Day" & group2 == "Stressed CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "1 Day", group2 == "Stressed CM"))



# Stressed 3 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "3 Days", Stress_Status == "Stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "3 Days" & group2 == "Stressed CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "3 Days", group2 == "Stressed CM"))



# Stressed 1 Week
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "1 Week", Stress_Status == "Stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "1 Week" & group2 == "Stressed CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "1 Week", group2 == "Stressed CM"))



# Stressed 3 Weeks
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "3 Weeks", Stress_Status == "Stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "3 Weeks" & group2 == "Stressed CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "3 Weeks", group2 == "Stressed CM"))





#############################################################################
# Not stressed
#############################################################################
# Not stressed 6 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "6 Hours", Stress_Status == "Not stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "6 Hours" & group2 == "Not stressed CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "6 Hours", group2 == "Not stressed CM"))



# Not stressed 12 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "12 Hours", Stress_Status == "Not stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "12 Hours" & group2 == "Not stressed CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "12 Hours", group2 == "Not stressed CM"))



# Not stressed 1 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "1 Day", Stress_Status == "Not stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "1 Day" & group2 == "Not stressed CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "1 Day", group2 == "Not stressed CM"))



# Not stressed 3 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "3 Days", Stress_Status == "Not stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "3 Days" & group2 == "Not stressed CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "3 Days", group2 == "Not stressed CM"))



# Not stressed 1 Week
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "1 Week", Stress_Status == "Not stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "1 Week" & group2 == "Not stressed CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "1 Week", group2 == "Not stressed CM"))



# Not Stressed 3 Weeks
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "3 Weeks", Stress_Status == "Not stressed CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "3 Weeks" & group2 == "Not stressed CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "3 Weeks", group2 == "Not stressed CM"))









#############################################################################
# SHAM
#############################################################################
# sham 6 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "6 Hours", Stress_Status == "SHAM CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "6 Hours" & group2 == "SHAM CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "6 Hours", group2 == "SHAM CM"))



# sham 12 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "12 Hours", Stress_Status == "SHAM CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "12 Hours" & group2 == "SHAM CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "12 Hours", group2 == "SHAM CM"))




# sham 1 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "1 Day", Stress_Status == "SHAM CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression"
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "1 Day" & group2 == "SHAM CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "1 Day", group2 == "SHAM CM"))



# sham 3 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "3 Days", Stress_Status == "SHAM CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "3 Days" & group2 == "SHAM CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "3 Days", group2 == "SHAM CM"))



# sham 1 Week
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "1 Week", Stress_Status == "SHAM CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "1 Week" & group2 == "SHAM CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "1 Week", group2 == "SHAM CM"))



# SHAM 3 Weeks
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
value <- go_process_summary %>%
  filter(gene == gene_of_interest, Timepoint == "3 Weeks", Stress_Status == "SHAM CM") %>%
  pull(mean_expression) 
print(value)

# Update "group2_mean_expression" 
anova_tukey_mean <- anova_tukey_mean %>%
  mutate(group2_mean_expression = ifelse(gene == gene_of_interest & Timepoint == "3 Weeks" & group2 == "SHAM CM", value, group2_mean_expression))

# Print to verify the update
print(anova_tukey_mean %>% filter(gene == gene_of_interest, Timepoint == "3 Weeks", group2 == "SHAM CM"))


view(anova_tukey_mean)
