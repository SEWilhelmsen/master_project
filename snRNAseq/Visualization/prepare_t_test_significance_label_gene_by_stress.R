# Prepare t_test_gene_by_stress for line plot
# Silje Wilhelmsen

# Load data
t_test_gene_by_stress <- read_excel("Data/Stress_Status/t_test_gene_by_stress.xlsx")
View(t_test_gene_by_stress)

# Add significance label if necessary in two columns
t_test_gene_by_stress <- t_test_gene_by_stress %>%
  mutate(significance_label = ifelse(group1 == "SHAM CM" & p_adj < 0.001, "***",
                                     ifelse(group1 == "SHAM CM" & p_adj < 0.01, "**",
                                            ifelse(group1 == "SHAM CM" & p_adj < 0.05, "*", "ns")))) %>%
  mutate(stressed_not_stressed = ifelse(group1 != "SHAM CM" & group2 != "SHAM CM" & p_adj < 0.05, "\u25B2", "ns"))
print(t_test_gene_by_stress)


gene_of_interest <- "IDH3B" # Change gene 

# Load gene summary data
# go_process_summary <- read_excel("Data/Stress_Status/go_process_summary_idh1.xlsx")
# go_process_summary <- read_excel("Data/Stress_Status/go_process_summary_idh2.xlsx")
# go_process_summary <- read_excel("Data/Stress_Status/go_process_summary_idh3b.xlsx")



# Filter the correct gene
filtered_t_test <- t_test_gene_by_stress %>%
  filter(gene == gene_of_interest) 
print(filtered_t_test)
print(go_process_summary)


# Remove labels from ANOVA
go_process_summary <- go_process_summary %>%
  dplyr::select(-significance_label, -stress_vs_not_stressed)

print(go_process_summary)

# Add labels from t_test
go_process_summary <- go_process_summary %>%
  left_join(filtered_t_test %>% dplyr::select(Timepoint, gene, group2, p_adj, significance_label, stressed_not_stressed), 
            by = c("Timepoint" = "Timepoint", "gene" = "gene", "Stress_Status" = "group2"))

print(go_process_summary)


# Save 
write.xlsx(go_process_summary, "Data/Stress_Status/go_process_summary_t_test_idh3b.xlsx") # Change gene name

