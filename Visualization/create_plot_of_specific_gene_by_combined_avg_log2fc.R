## Create trend graph from combined_avg_log2fc (*ikke ferdig*)
# Input:combined_avg_log2fc data frame
# Output: plot of gene expression over time
# Silje Wilhelmsen

library(ggplot2)
library(dplyr)
library(tidyr)

### Load list of specific genes 
source("find_specific_genes.R")


# Filter for "ACTA1"
acta1_data <- common_genes_with_avg_log2fc %>%
  filter(gene == "ACTA1")

# Create a mapping of timepoint labels 
acta1_data <- acta1_data %>%
  rename("1 day" = "avg_log2fc_24h",
         "3 days" = "avg_log2fc_3d",
         "1 week" = "avg_log2fc_1w",
         "3 weeks" = "avg_log2fc_3w"
         )

head(acta1_data)

# Reshape data from wide to long format
long_format_acta1 <- acta1_data %>%
  pivot_longer(cols = -gene, 
               names_to = "Timepoint", 
               values_to = "Expression")

# Print to confirm the structure
print(head(long_format_acta1))

# Ensure the Timepoint factor levels are in the desired order
long_format_acta1$Timepoint <- factor(long_format_acta1$Timepoint,
                                      levels = c("1 day",
                                                 "3 days",
                                                 "1 week",
                                                 "3 weeks"))

# # Create a mapping of timepoint labels to days
# timepoint_mapping <- c(
#   "avg_log2FC_24h" = 1,   # 1 day
#   "avg_log2FC_3d" = 3,    # 3 days
#   "avg_log2FC_1w" = 7,    # 1 week (7 days)
#   "avg_log2FC_3w" = 21    # 3 weeks (21 days)
# )


# # Add numerical time values based on the mapping
# long_format_acta1 <- long_format_acta1 %>%
#   mutate(Day = timepoint_mapping[Timepoint])

# Check if there are any missing values
print(any(is.na(long_format_acta1)))

# # Ensure the Timepoint factor levels are in the desired order
# long_format_acta1$Timepoint <- factor(long_format_acta1$Timepoint, levels = names(Timepoint))


# Create the ggplot for ACTA1 with conditional coloring, numerical x-axis, and axis lines forming a frame
ggplot(long_format_acta1, aes(x = Timepoint, y = Expression, group = gene, color = gene)) +
  geom_line(size = 2) +
  geom_point(size = 2) +
  scale_color_manual(values = c("ACTA1" = "red")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  xlab("Time Point") +
  ylab("Average log2 Fold Change Expression") +
  ggtitle("Average log2FC of Expression of ACTA1 over Time") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = 'black')
  )



