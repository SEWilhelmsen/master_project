# Create plot of celltype proportion relative to the mean proportion value of that celltype or
# the mean of all celltypes


library(ggplot2)
library(dplyr)
library(cowplot) 

cell_name <- "Mesothelial"
meso_weights <- norm_weights[, cell_name]

# Estimate mean
# Alternative A: mean of all spots
mean_all <- mean(meso_weights, na.rm = TRUE)
# Alternative B: mean of spots with a value >0 
mean_nonzero <- mean(meso_weights[meso_weights > 0], na.rm = TRUE)

# Set mean to use
mean_val <- mean_all   # mean_nonzero or mean_all

# Create gradient-data from 0 to mean
meso_for_grad <- pmin(meso_weights, mean_val)

# Create plot of background gradient
test_plot <- try(plot_puck_continuous(RCTD_temp@spatialRNA, barcodes, meso_for_grad,
                                      title = paste0(cell_name, " (gradient towards mean, >mean = red) ", "(", angle_label, ")"),
                                      size = 1), silent = TRUE)

if (inherits(test_plot, "ggplot")) {
  base_plot <- test_plot + coord_fixed(xlim = xlim_all, ylim = ylim_all, expand = FALSE)
} else {
  plot_puck_continuous(RCTD_temp@spatialRNA, barcodes, meso_for_grad,
                       title = paste0(cell_name, " (gradient towards mean, >mean = red) ", "(", angle_label, ")"),
                       size = 1,
                       xlim = xlim_all, ylim = ylim_all)
  base_plot <- NULL
}

# Create another dataset with spots > mean
sp <- RCTD_temp@spatialRNA
coords_df <- if (!is.null(sp@coords)) as.data.frame(sp@coords) else as.data.frame(sp@coordinates)
colnames(coords_df)[1:2] <- c("x", "y")

df <- data.frame(barcode = barcodes, meso = meso_weights, stringsAsFactors = FALSE)
coords_matched <- coords_df[match(df$barcode, rownames(coords_df)), c("x","y")]
df_coords <- bind_cols(df, coords_matched)

df_red <- df_coords %>% filter(!is.na(meso) & meso > mean_val)

# Overlay red points
if (!is.null(base_plot)) {
  p <- base_plot +
    geom_point(data = df_coords %>% filter(!is.na(x) & !is.na(y)),
               aes(x = x, y = y),
               color = "grey90", size = 0.6) +
    geom_point(data = df_red, aes(x = x, y = y), color = "red", size = 2, alpha = 0.9) +
    ggtitle(paste0(cell_name, " (mean = ", signif(mean_val, 3), ", sample: ", prefix, ")"))
  print(p)
} else {
  points(df_red$x, df_red$y, col = "red", pch = 16, cex = 1.2)
}
