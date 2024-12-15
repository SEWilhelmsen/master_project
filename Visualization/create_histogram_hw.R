# Create histogram 
# Input: table (excel?) of animal data
# Output: histograms of all samples, by condition and time_point
# Silje Wilhelmsen

# Create histograms separated by condition and time_point
histogram_hw <- processed_data %>%
  ggplot(aes(x = `hw`, fill = condition)) +
  geom_histogram(binwidth = 6, color = "black", alpha = 0.7, position = "dodge") +
  facet_grid(condition ~ time_point) +
  labs(title = "Heart Weight of all samples", 
       x = "HW (mg)", 
       y = "Frequency") +
  scale_fill_manual(values = c("sham" = "white", "AB" = "coral")) + 
  theme_minimal() +
  theme(legend.position = "none")

# Display the plot
print(histogram_hw)

# # Save
# png(filename = "Y:/Silje/snRNAseq/Plots/histo_all_samples_hw.png", width = 1000, height = 800)
