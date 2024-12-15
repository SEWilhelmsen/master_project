## Create a grouped barplot from a data frame
# Silje Wilhelmsen


library(ggplot2)

# Create data - How to import the real data??? 
data <- data.frame(
  timepoint = rep(c("24 hours", "3 days", "1 week", "3 weeks"), each = 2), 
  condition = rep(c("AB", "SHAM"), times = 4), 
  value = c(5, 3, 6, 4, 7, 5, 2, 6)
)

data$timepoint <- factor(data$timepoint, 
                         levels = c("24 hours", "3 days", "1 week", "3 weeks"))

ggplot(data, aes(x = timepoint, y = value, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Bar plot of conditions over timepoints", 
       x = "Timepoint",
       y = "Value") + 
  theme_minimal() +
  scale_fill_manual(values = c("AB" = "blue", "SHAM" = "coral"))