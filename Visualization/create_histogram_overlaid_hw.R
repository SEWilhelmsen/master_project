# Create histogram with overlaid graphs
# Input: table (excel?) of animal data
# Output: histograms of all samples, by condition and time_point
# Silje Wilhelmsen

# Create overlaid histograms separated by condition and time_point
histogram_overlay_hw <- animal_data_rnaseq_renamed %>%
  ggplot(aes(x = `HW_(mg)`, fill = condition, color = condition)) +
  geom_histogram(binwidth = 6, alpha = 0.4, position = "identity") +
  facet_grid(. ~ time_point) +
  labs(title = "Heart Weight of All Samples", 
       x = "HW (mg)", 
       y = "Frequency") +
  scale_fill_manual(values = c("sham" = "blue", "AB" = "coral")) + 
  scale_color_manual(values = c("sham" = "blue", "AB" = "coral")) +
  theme_minimal() +
  theme(legend.position = "right")


# Display the plot
print(histogram_overlay_hw)

# Save the plot 
png(filename = "Y:/Silje/snRNAseq/Plots/histo_all_samples_hw.png", width = 1000, height = 800)
