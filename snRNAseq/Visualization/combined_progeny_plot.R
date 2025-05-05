# Create one plot consisting of six PROGENy figures


# Make a figure consisting of multiple figures
##############################################################################################

progeny_6h <- progeny_6h +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin
progeny_12h <- progeny_12h +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin

progeny_1d <- progeny_1d +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin
progeny_3d <- progeny_3d +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin

progeny_1w <- progeny_1w +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin
progeny_3w <- progeny_3w +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin


# Example of laying out in a 3x3 grid
total_progeny_plot <- (progeny_6h | progeny_12h) /
  (progeny_1d | progeny_3d) /
  (progeny_1w | progeny_3w)

# Print the combined plot
print(total_progeny_plot)

output_dir_progeny <- "C:/Users/siljeew/Master_project/snRNAseq/Plots/PROGENy/"
# Save the combined plot—adjust width and height for space between plots
ggsave(file.path(output_dir_progeny, "total_progeny_plot.pdf"), plot = total_progeny_plot, width = 20, height = 22)

ggsave(file.path(output_dir_progeny, "total_progeny_plot.png"), plot = total_progeny_plot, width = 20, height = 22)

  