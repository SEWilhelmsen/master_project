# Split spatial transcriptomics slide into four aectors
# Forutsetninger:
# - 'object' er ditt Seurat-objekt
# - evt. 'RCTD_temp' objekt finnes (hvis du brukte RCTD)
# Use these guidelines for segmentation: https://www.ahajournals.org/doi/10.1161/hc0402.102975


library(dplyr)
library(ggplot2)

# -------------------------------------------------------------------------
# Get coordinates
# -------------------------------------------------------------------------
coords <- NULL
if (exists("RCTD") && !is.null(RCTD) && !is.null(RCTD@spatialRNA@coords)) {
  coords_df <- as.data.frame(RCTD@spatialRNA@coords)
  coords_df$spot <- rownames(coords_df)
}

# Standardiser kolonnenavn til x,y (tilpass hvis kolonner heter annet)
colnames(coords_df)[1:2] <- c("x", "y")
coords <- coords_df
head(coords)
head(coords_df)

# -------------------------------------------------------------------------
# Define centre
# -------------------------------------------------------------------------
# Set centre manually
xc <- 4500; yc <- 4500 # TP48.1
xc <- 4500; yc <- 5500 # TP36.1
xc <- 5600; yc <- 5000 # TP37.5
xc <- 5600; yc <- 4200 # TPA1

# Angi vinkler i grader som du vil vise 
angles_deg <- c(0, 60, 120, 200, 260, 320)   # TP36.1
angles_deg <- c(10, 70, 130, 200, 270, 320)   # TP37.5
angles_deg <- c(0, 70, 140, 200, 260, 310)   # TPA1

radius <- max(diff(range(coords$x)), diff(range(coords$y))) * 0.8

# Lag linjedata for hver vinkel
rays <- do.call(rbind, lapply(angles_deg, function(a) {
  rad <- a * pi / 180
  xend <- xc + radius * cos(rad)
  yend <- yc + radius * sin(rad)
  data.frame(x = xc, y = yc, xend = xend, yend = yend, angle = a)
}))

ggplot(coords, aes(x = x, y = y)) +
  geom_point(size = 0.5, color = "grey60") +
  geom_segment(data = rays, aes(x = x, y = y, xend = xend, yend = yend),
               color = "blue", arrow = arrow(length = unit(0.02, "npc")), linewidth = 0.8) +
  geom_point(aes(x = xc, y = yc), color = "red", size = 3) +
  coord_equal() + theme_minimal() +
  ggtitle("Donde esta el centro?")



# # Set centre based on mean
# center_x <- mean(coords$x)
# center_y <- mean(coords$y)
# xc <- center_x; yc <- center_y
# 
# # Angi vinkler i grader som du vil vise (eks. kvartaler)
# angles_deg <- c(0, 90, 180, 270)   # du kan endre til dine sektor-grenser, f.eks. c(330,30,120,210)
# radius <- max(diff(range(coords$x)), diff(range(coords$y))) * 0.8
# 
# # Lag linjedata for hver vinkel
# rays <- do.call(rbind, lapply(angles_deg, function(a) {
#   rad <- a * pi / 180
#   xend <- xc + radius * cos(rad)
#   yend <- yc + radius * sin(rad)
#   data.frame(x = xc, y = yc, xend = xend, yend = yend, angle = a)
# }))
# 
# ggplot(coords, aes(x = x, y = y)) +
#   geom_point(size = 0.5, color = "grey60") +
#   geom_segment(data = rays, aes(x = x, y = y, xend = xend, yend = yend),
#                color = "blue", arrow = arrow(length = unit(0.02, "npc")), size = 0.8) +
#   geom_point(aes(x = xc, y = yc), color = "red", size = 3) +
#   coord_equal() + theme_minimal() +
#   ggtitle("Donde esta el centro?")

# -------------------------------------------------------------------------
# Estimate angles
# -------------------------------------------------------------------------
coords <- coords %>%
  mutate(dx = x - xc,
         dy = y - yc,
         theta_rad = atan2(dy, dx),
         theta_deg = (theta_rad * 180 / pi) %% 360)  # 0..360

# Hvis du trenger å rotere hele systemet (offset) for anatomisk orientering:
offset_deg <- 0    # sett f.eks. 90 eller 45 hvis du trenger rotasjon
coords$theta_deg_rot <- (coords$theta_deg + offset_deg) %% 360

# -------------------------------------------------------------------------
# Define sectors
# -------------------------------------------------------------------------
# TP36.1
# sector_defs <- data.frame(
#   region = c("Mid anterior", "Mid anterior septum", "Septum", "Mid inferior", "Mid post", "Mid lateral"),
#   start  = c(0, 60, 120, 200, 260, 320) ,   # eksempel; se forklaring under
#   end    = c(60, 120, 200, 260, 320, 0)    # exclusive end
# )


# # TP37.5
# sector_defs <- data.frame(
#   region = c("Mid anterior", "Mid anterior septum", "Septum", "Mid inferior", "Mid post", "Mid lateral"),
#   start  = c(10, 70, 130, 200, 270, 320) ,   # eksempel; se forklaring under
#   end    = c(70, 130, 200, 270, 320, 10)    # exclusive end
# )

# A1
sector_defs <- data.frame(
  region = c("Mid anterior", "Mid anterior septum", "Septum", "Mid inferior", "Mid post", "Mid lateral"),
  start  = c(0, 70, 140, 200, 260, 310) ,   # eksempel; se forklaring under
  end    = c(70, 140, 200, 260, 310, 0)    # exclusive end
)

# NB: rekkefølgen i region/start/end må samsvare. Du kan bytte grenser slik du ønsker.

angle_to_region <- function(angle_deg, defs) {
  res <- rep(NA_character_, length(angle_deg))
  for (i in seq_len(nrow(defs))) {
    s <- defs$start[i]; e <- defs$end[i]; name <- defs$region[i]
    if (s < e) {
      inside <- angle_deg >= s & angle_deg < e
    } else { # wrap-around
      inside <- (angle_deg >= s & angle_deg < 360) | (angle_deg >= 0 & angle_deg < e)
    }
    res[inside] <- name
  }
  return(res)
}

# Bruk rotert vinkel (theta_deg_rot)
coords$Sector <- angle_to_region(coords$theta_deg_rot, sector_defs)


# -------------------------------------------------------------------------
# Visualize
# -------------------------------------------------------------------------
print(table(coords$Sector, useNA = "ifany"))

coords$Sector <- factor(coords$Sector, levels = c("Mid anterior", "Mid anterior septum", "Septum",  "Mid inferior", "Mid post", "Mid lateral"))

ggplot(coords, aes(x = x, y = y, color = Sector)) +
  geom_point(size = 0.8) + coord_equal() +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  ggtitle(paste0("Sectors (sample: ", prefix, ")")) +
  theme_minimal()

# -------------------------------------------------------------------------
# Add sector information to meta.data 
# -------------------------------------------------------------------------
df_sectors <- coords %>% select(spot, Sector)
match_idx <- match(rownames(object@meta.data), df_sectors$spot)
object@meta.data$Sector <- df_sectors$Sector[match_idx]
desired_levels <- c("Septum", "Free wall", "Hinge 1", "Hinge 2")
# Hvis noen labels i df_sectors ikke matcher ønsket_levels, remap før faktor hvis nødvendig.
# object@meta.data$Sector <- factor(object@meta.data$Sector, levels = desired_levels)
print(table(object@meta.data$Sector, useNA = "ifany"))

# -------------------------------------------------------------------------
# Quality control
# -------------------------------------------------------------------------
head(rownames(object@meta.data))
head(df_sectors$spot)
length(intersect(rownames(object@meta.data), df_sectors$spot))
head(coords)  



# -------------------------------------------------------------------------
# Define concentric zones
# -------------------------------------------------------------------------

if (!"spot" %in% colnames(coords)) coords$spot <- rownames(coords)

# 1) dist
coords <- coords %>% mutate(dist = sqrt((x - xc)^2 + (y - yc)^2))

# 2) thresholds (bruk faste som i ditt eksempel)
n_rings <- 3
max_dist <- max(coords$dist, na.rm = TRUE)
thresholds <- c(0, 0.2 * max_dist, 0.72 * max_dist, 1.00 * max_dist)


# guard: sørg for sorterte og unike terskler
thresholds <- sort(unique(thresholds))
if (length(thresholds) != (n_rings + 1)) {
  warning("Tersklene er ikke av forventet lengde (mulig duplikater).")
}

# 3) zone names
zone_names <- c("Central", "Mid", "Outer")
if (length(zone_names) != n_rings) stop("Lengde på zone_names må være lik n_rings")

# 4) assign_zone funksjon (som du har)
assign_zone <- function(d, thr, names) {
  k <- length(thr) - 1
  res <- rep(NA_character_, length(d))
  for (i in seq_len(k)) {
    lo <- thr[i]; hi <- thr[i+1]
    inside <- d >= lo & d <= hi
    if (i < k) inside <- d >= lo & d < hi
    res[inside] <- names[i]
  }
  res
}

# 5) assign
coords$Zone <- assign_zone(coords$dist, thresholds, zone_names)
coords$Zone[is.na(coords$Zone) & coords$dist <= max(coords$dist, na.rm = TRUE)] <- zone_names[length(zone_names)]

# 6) push to Seurat
df_zones <- coords %>% select(spot, Zone)
match_idx <- match(rownames(object@meta.data), df_zones$spot)
object@meta.data$Zone <- df_zones$Zone[match_idx]
object@meta.data$Zone <- factor(object@meta.data$Zone, levels = zone_names)

# 7) checks
print(table(coords$Zone, useNA = "ifany"))
print(sum(!is.na(match_idx)))

# 8) plot (som du har)
theta <- seq(0, 2*pi, length.out = 360)
circle_df <- do.call(rbind, lapply(thresholds[-1], function(r) {
  data.frame(x = xc + r * cos(theta), y = yc + r * sin(theta), radius = r)
}))

ggplot(coords, aes(x = x, y = y, color = Zone)) +
  geom_point(size = 1.5, alpha = 0.9) +
  geom_path(data = circle_df, aes(x = x, y = y, group = radius),
            color = "black", size = 0.6, linetype = "dashed") +
  geom_point(aes(x = xc, y = yc), color = "red", size = 3) +
  coord_equal() +
  scale_color_brewer(palette = "Spectral") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  ggtitle("Concentric zones")




# 9) (Valgfritt) Kombiner sektor + sone til en felles label
coords <- coords %>% mutate(SectorZone = paste0(Sector, "_", Zone))
df_sz <- coords %>% select(spot, SectorZone)
match_idx2 <- match(rownames(object@meta.data), df_sz$spot)
object@meta.data$SectorZone <- df_sz$SectorZone[match_idx2]

# Siste sjekk:
print(table(object@meta.data$SectorZone, useNA = "ifany"))


theta <- seq(0, 2*pi, length.out = 360)
circle_df <- do.call(rbind, lapply(thresholds[-1], function(r) {
  data.frame(x = xc + r * cos(theta),
             y = yc + r * sin(theta),
             radius = r)
}))

# Lag sector-ray data (valgfritt). Hvis du har sektor-grenser i grader (f.eks. sector_defs$start),
# tegn rays ved start-vinklene:
# Eksempel: sector_starts <- c(330, 30, 120, 210)
sector_starts <- if (exists("sector_defs")) sector_defs$start else c(0,90,180,270)
radius_max <- max(thresholds, na.rm = TRUE) * 1.05
rays <- do.call(rbind, lapply(sector_starts, function(a) {
  rad <- a * pi / 180
  data.frame(x = xc, y = yc,
             xend = xc + radius_max * cos(rad),
             yend = yc + radius_max * sin(rad),
             angle = a)
}))

# --------------------------------------------------------------------------
# A) Farge = Sector, tegn ringene (best for å se sektorer med sone-grenser)
# --------------------------------------------------------------------------
pA <- ggplot(coords, aes(x = x, y = y, color = Sector)) +
  geom_point(size = 1.4, alpha = 0.9) +
  geom_path(data = circle_df, aes(x = x, y = y, group = radius),
            color = "black", size = 0.6, linetype = "dashed") +
  geom_segment(data = rays, aes(x = x, y = y, xend = xend, yend = yend),
               color = "grey30", size = 0.6) +
  geom_point(aes(x = xc, y = yc), color = "red", size = 3) +
  coord_equal() +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  ggtitle("Sectors (color) + concentric zones (dashed circles)")

print(pA)



# --------------------------------------------------------------------------
# B) Facet per Zone (viser sektorfordeling innen hver sone)
# --------------------------------------------------------------------------
pB <- ggplot(coords, aes(x = x, y = y, color = Sector)) +
  geom_point(size = 1.2, alpha = 0.9) +
  geom_path(data = circle_df, aes(x = x, y = y, group = radius),
            color = "black", size = 0.5, linetype = "dashed") +
  geom_point(aes(x = xc, y = yc), color = "red", size = 2) +
  coord_equal() +
  facet_wrap(~ Zone) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  ggtitle("Facet per Zone: Sectors within each concentric ring")

print(pB)


# --------------------------------------------------------------------------
# C) Color = SectorZone (hver kombinasjon får egen farge) — kan bli mange farger
# --------------------------------------------------------------------------
pC <- ggplot(coords, aes(x = x, y = y, color = SectorZone)) +
  geom_point(size = 1.6, alpha = 0.95) +
  geom_path(data = circle_df, aes(x = x, y = y, group = radius),
            color = "black", size = 0.5, linetype = "dashed") +
  coord_equal() +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 4))) +
  theme_minimal() +
  ggtitle("SectorZone (combined)")

print(pC)


