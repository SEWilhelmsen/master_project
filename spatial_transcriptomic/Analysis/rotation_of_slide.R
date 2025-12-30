# Rotate slide and save if necessary as part of run_create_plots_of_cell_type_location.R
# Hardkodet rotasjon i grader (positiv = CCW). Sett ønsket vinkel her:
angle_rotation <- -60   #  -90 for 90 degrees clockwise, 90 degrees for counter clockwise
angle_rotation <- 0

# Hjelpefunksjon for generell rotasjon (returnerer data.frame, bevarer navn)
rotate_any <- function(coords, angle_deg, center = NULL) {
  coords <- as.data.frame(coords)
  if (is.null(center)) center <- c(mean(coords[[1]]), mean(coords[[2]]))
  theta <- angle_deg * pi / 180
  R <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2, byrow = TRUE)
  mat <- as.matrix(coords)
  shifted <- sweep(mat, 2, center, "-")
  rotated <- t(R %*% t(shifted))
  rotated <- sweep(rotated, 2, center, "+")
  res <- as.data.frame(rotated)
  colnames(res) <- colnames(coords)
  rownames(res) <- rownames(coords)
  res
}

# Hent original coords som data.frame
coords_orig_df <- if (is.data.frame(coords)) coords else as.data.frame(coords)
center <- c(mean(coords_orig_df[[1]]), mean(coords_orig_df[[2]]))

# Roter (hvis angle_rotation == 0, coords_rot_df == coords_orig_df)
if (angle_rotation == 0) {
  coords_rot_df <- coords_orig_df
  angle_label <- "none"
} else {
  coords_rot_df <- rotate_any(coords_orig_df, angle_rotation, center = center)
  angle_label <- paste0("rot", angle_rotation)
}

# Lag en kopi av RCTD og sett inn roterte coords i spatialRNA (beholder originalen)
RCTD_temp <- RCTD
sp_temp <- RCTD_temp@spatialRNA

if (!is.null(sp_temp@coords)) {
  sp_temp@coords <- coords_rot_df
} else if (!is.null(sp_temp@coordinates)) {
  sp_temp@coordinates <- coords_rot_df
} else {
  stop("Fant verken @coords eller @coordinates i spatialRNA-objektet")
}
RCTD_temp@spatialRNA <- sp_temp

# Preview
i <- 19
cell_preview <- cell_type_names[i]

test_plot <- try(plot_puck_continuous(RCTD_temp@spatialRNA, barcodes, norm_weights[, cell_preview],
                                      title = paste0(cell_preview, " (", angle_label, ")"), size = 1), silent = TRUE)

if (inherits(test_plot, "ggplot")) {
  print(test_plot + coord_fixed(xlim = xlim_all, ylim = ylim_all, expand = FALSE))
} else {
  plot_puck_continuous(RCTD_temp@spatialRNA, barcodes, norm_weights[, cell_preview],
                       title = paste0(cell_preview, " (", angle_label, ")"), size = 1,
                       xlim = xlim_all, ylim = ylim_all)
}

# Save RCTD_temp
data_dir <- "C:/Users/Labuser/master_project/Spatials/tmp/"
saveRDS(RCTD_temp, file = paste0(data_dir, prefix, "RCTD_temp.rds"))


# Dersom du er fornøyd med rotasjonen kan du lagre dette for alle celletypene
test_plot2 <- try(plot_puck_continuous(RCTD_temp@spatialRNA, barcodes, norm_weights[, cell_type_names[1]],
                                       title = cell_type_names[1], size = 1), silent = TRUE)
ggplot_mode <- inherits(test_plot2, "ggplot")

if (!exists("prefix")) prefix <- ""



for (j in seq_len(n_types)) {
  cellname <- cell_type_names[j]
  fname <- file.path(plot_dir, paste0(prefix, "distribution_", j, "_", cellname, angle_rotation, ".tiff"))
  if (ggplot_mode) {
    p <- plot_puck_continuous(RCTD_temp@spatialRNA, barcodes, norm_weights[, cellname], title = cellname, size = 1)
    p_fixed <- p + coord_fixed(xlim = xlim_all, ylim = ylim_all, expand = FALSE)
    ggsave(filename = fname, plot = p_fixed, device = "tiff", width = 10, height = 8, units = "in", dpi = 600)
  } else {
    tiff(filename = fname, width = 10, height = 8, units = "in", res = 600)
    op <- par(pty = "s", mar = c(3,3,3,3))
    plot_puck_continuous(RCTD_temp@spatialRNA, barcodes, norm_weights[, cellname], title = cellname, size = 1,
                         xlim = xlim_all, ylim = ylim_all, asp = 1)
    par(op); dev.off()
  }
  message("Lagret: ", fname)
}

message("Ferdig. Rotasjon brukt: ", angle_label)

