# Create histogram (ikke ferdig)
# Input: table (excel?) of animal data
# Output: histograms of all samples, by condition and time_point
# Silje Wilhelmsen

# # Load libraries
# library(ggplot2)




# Create histogram
histogram_variable <- processed_data %>%
  ggplot(processed_data, aes(x = variable_of_interest)) +
  geom_histogram(binwidth = 1, color = "white")


# Display the plot
print(histogram_variable)

# Save
png(filename = "Y:/Silje/snRNAseq/Plots/histo_all_samples_hw.png", width = 1000, height = 800)
