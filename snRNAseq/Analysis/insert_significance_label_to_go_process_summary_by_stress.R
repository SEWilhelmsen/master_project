# Insert significance labels to go_process_summary

# go_process_summary <- go_process_summary %>%
#   dplyr::select(-significance_label, -significance_label_group2.x, -significance_label_group2.y, -gene.y, -gene.y.y, -gene.x.x, -gene.x)


# Create column for significance labels
go_process_summary <- go_process_summary %>%
  mutate(significance_label = NA)

# Pull the value from anova_tukey_means and insert into go_process_summary
# Stressed 6 Hours compared to SHAM
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "6 Hours", group1 == "Stressed CM", group2 == "SHAM CM") %>%
  pull(significance_label) 

# Print the retrieved value
print(padj_label)


# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(gene == gene_of_interest & Timepoint == "6 Hours" & Stress_Status == "Stressed CM", padj_label, significance_label))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "6 Hours", Stress_Status == "Stressed CM"))



# Stressed 12 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "12 Hours", group1 == "Stressed CM", group2 == "SHAM CM") %>%
  pull(significance_label) 
print(padj_label)

# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(gene == gene_of_interest & Timepoint == "12 Hours" & Stress_Status == "Stressed CM", padj_label, significance_label))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "12 Hours", Stress_Status == "Stressed CM"))



# Stressed 1 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "1 Day", group1 == "Stressed CM", group2 == "SHAM CM") %>%
  pull(significance_label) 

# Print the retrieved value
print(padj_label)

# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(gene == gene_of_interest & Timepoint == "1 Day" & Stress_Status == "Stressed CM", padj_label, significance_label))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "1 Day", Stress_Status == "Stressed CM"))



# Stressed 3 Days
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "3 Days", group1 == "Stressed CM", group2 == "SHAM CM") %>%
  pull(significance_label) 
print(padj_label)

# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(gene == gene_of_interest & Timepoint == "3 Days" & Stress_Status == "Stressed CM", padj_label, significance_label))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "3 Days", Stress_Status == "Stressed CM"))



# Stressed 1 Week
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "1 Week", group1 == "Stressed CM", group2 == "SHAM CM") %>%
  pull(significance_label) 
print(padj_label)

# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(gene == gene_of_interest & Timepoint == "1 Week" & Stress_Status == "Stressed CM", padj_label, significance_label))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "1 Week", Stress_Status == "Stressed CM"))



# Stressed 3 Weeks
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "3 Weeks", group1 == "Stressed CM", group2 == "SHAM CM") %>%
  pull(significance_label) 
print(padj_label)

# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(gene == gene_of_interest & Timepoint == "3 Weeks" & Stress_Status == "Stressed CM", padj_label, significance_label))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "3 Weeks", Stress_Status == "Stressed CM"))






# Not stressed 6 Hours compared to SHAM
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "6 Hours", group2 == "Not stressed CM", group1 == "SHAM CM") %>%
  pull(significance_label) 
print(padj_label)

# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(gene == gene_of_interest & Timepoint == "6 Hours" & Stress_Status == "Not stressed CM", padj_label, significance_label))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "6 Hours", Stress_Status == "Not stressed CM"))



# Stressed 12 Hours
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "12 Hours", group2 == "Not stressed CM", group1 == "SHAM CM") %>%
  pull(significance_label) 
print(padj_label)


# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(gene == gene_of_interest & Timepoint == "12 Hours" & Stress_Status == "Not stressed CM", padj_label, significance_label))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "12 Hours", Stress_Status == "Not stressed CM"))



# Stressed 1 Day
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "1 Day", group2 == "Not stressed CM", group1 == "SHAM CM") %>%
  pull(significance_label) 
print(padj_label)

# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(gene == gene_of_interest & Timepoint == "1 Day" & Stress_Status == "Not stressed CM", padj_label, significance_label))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "1 Day", Stress_Status == "Not stressed CM"))



# Stressed 3 Days
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "3 Days", group2 == "Not stressed CM", group1 == "SHAM CM") %>%
  pull(significance_label) 

# Print the retrieved value
print(padj_label)

# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(gene == gene_of_interest & Timepoint == "3 Days" & Stress_Status == "Not stressed CM", padj_label, significance_label))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "3 Days", Stress_Status == "Not stressed CM"))



# Stressed 1 Week
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "1 Week", group1 == "SHAM CM", group2 == "Not stressed CM", group1 == "SHAM CM") %>%
  pull(significance_label) 
print(padj_label)

# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(gene == gene_of_interest & Timepoint == "1 Week" & Stress_Status == "Not stressed CM", padj_label, significance_label))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "1 Week", Stress_Status == "Not stressed CM"))



# Stressed 3 Weeks
#############################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "3 Weeks", group2 == "Not stressed CM", group1 == "SHAM CM") %>%
  pull(significance_label) 
print(padj_label)

# Update significance_label
go_process_summary <- go_process_summary %>%
  mutate(significance_label = ifelse(gene == gene_of_interest & Timepoint == "3 Weeks" & Stress_Status == "Not stressed CM", padj_label, significance_label))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "3 Weeks", Stress_Status == "Not stressed CM"))



View(go_process_summary)









####################################################################################
# Create new column for labeling if significant between stressed and not stressed
###################################################################################
go_process_summary <- go_process_summary %>%
  mutate(stress_vs_not_stressed = NA)


# Stressed 6 Hours compared to not stressed
####################################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "6 Hours", group1 == "Stressed CM", group2 == "Not stressed CM") %>%
  pull(significance_label) 
print(padj_label)


# Update stress_vs_not_stressed
go_process_summary <- go_process_summary %>%
  mutate(stress_vs_not_stressed = ifelse(gene == gene_of_interest & Timepoint == "6 Hours" & Stress_Status == "Stressed CM", padj_label, stress_vs_not_stressed))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "6 Hours", Stress_Status == "Stressed CM"))



# Stressed 12 Hours
###################################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "12 Hours", group1 == "Stressed CM", group2 == "Not stressed CM") %>%
  pull(significance_label) 
print(padj_label)

# Update stress_vs_not_stressed
go_process_summary <- go_process_summary %>%
  mutate(stress_vs_not_stressed = ifelse(gene == gene_of_interest & Timepoint == "12 Hours" & Stress_Status == "Stressed CM", padj_label, stress_vs_not_stressed))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "12 Hours", Stress_Status == "Stressed CM"))



# Stressed 1 Day
##################################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "1 Day", group1 == "Stressed CM", group2 == "Not stressed CM") %>%
  pull(significance_label) 
print(padj_label)

# Update stress_vs_not_stressed
go_process_summary <- go_process_summary %>%
  mutate(stress_vs_not_stressed = ifelse(gene == gene_of_interest & Timepoint == "1 Day" & Stress_Status == "Stressed CM", padj_label, stress_vs_not_stressed))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "1 Day", Stress_Status == "Stressed CM"))



# Stressed 3 Days
##################################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "3 Days", group1 == "Stressed CM", group2 == "Not stressed CM") %>%
  pull(significance_label) 
print(padj_label)

# Update stress_vs_not_stressed
go_process_summary <- go_process_summary %>%
  mutate(stress_vs_not_stressed = ifelse(gene == gene_of_interest & Timepoint == "3 Days" & Stress_Status == "Stressed CM", padj_label, stress_vs_not_stressed))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "3 Days", Stress_Status == "Stressed CM"))



# Stressed 1 Week
##################################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "1 Week", group1 == "Stressed CM", group2 == "Not stressed CM") %>%
  pull(significance_label) 
print(padj_label)

# Update stress_vs_not_stressed
go_process_summary <- go_process_summary %>%
  mutate(stress_vs_not_stressed = ifelse(gene == gene_of_interest & Timepoint == "1 Week" & Stress_Status == "Stressed CM", padj_label, stress_vs_not_stressed))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "1 Week", Stress_Status == "Stressed CM"))



# Stressed 3 Weeks
##################################################################################
# Filter to get specific row and print the mean expression for the given conditions
padj_label <- anova_tukey_mean %>%
  filter(gene == gene_of_interest, Timepoint == "3 Weeks", group1 == "Stressed CM", group2 == "Not stressed CM") %>%
  pull(significance_label) 
print(padj_label)

# Update stress_vs_not_stressed
go_process_summary <- go_process_summary %>%
  mutate(stress_vs_not_stressed = ifelse(gene == gene_of_interest & Timepoint == "3 Weeks" & Stress_Status == "Stressed CM", padj_label, stress_vs_not_stressed))

# Print to verify the update
print(go_process_summary %>% filter(gene == gene_of_interest, Timepoint == "3 Weeks", Stress_Status == "Stressed CM"))



View(go_process_summary)

