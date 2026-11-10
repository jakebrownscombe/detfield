# Title: Quadrature error of location_surface() against its node count
# Description: The numbers that ?location_surface states under Computation
#   about the node count. On the shipped example (stoney_track, 70
#   receivers, 4,818 water cells) with a receiver-hour spread of 2.55
#   logits, surfaces from 1, 20 and 64 Gauss-Hermite nodes are compared with
#   surfaces from 200 nodes, by total variation per animal-hour. The time
#   per animal is in validation/03_location_surface_timing.R.
# Run from the package root: Rscript validation/01_node_error.R
# Output: validation/outputs/01_node_error.csv and the printed lines below,
#   kept in validation/outputs/01_node_error.txt.

# ---- Libraries ----
library(detfield)

# ---- Field and counts of the example ----
lake <- terra::rast(system.file("extdata", "stoney_depth.tif", package = "detfield"))
truth <- stoney_track$truth
curve <- de_logistic(d50 = truth$d50, d05 = truth$d05, depth = truth$depth)
field <- de_field(water_distances(lake, stoney_track$receivers),
                  de_model(curve, sigma = list(station_hour = 2.55)))
counts <- bin_detections(stoney_track$detections, stoney_track$receivers, stoney_track$tags)
length(field$cells)  # water cells
nrow(unique(counts[c("animal_id", "bin_start")]))  # animal-hours with a detection

# ---- Surfaces at 200 nodes, the reference ----
# kept_mass = 1 keeps every water cell, so the comparison covers whole surfaces.
full <- function(nodes) {
  suppressWarnings(location_surface(counts, field,
                                    control = list(nodes = nodes, kept_mass = 1)))
}
ref <- full(200)
key <- function(s) paste(s$animal_id, format(s$bin_start, "%Y-%m-%d %H:%M"), s$cell)
length(unique(paste(ref$animal_id, ref$bin_start)))  # animal-hours in the reference

# ---- Total variation per animal-hour ----
# Half the summed absolute difference of the two surfaces over every water
# cell. A cell kept in one surface and not the other counts with 0 in the other.
tv_per_hour <- function(s) {
  keys <- union(key(s), key(ref))
  a <- s$prob[match(keys, key(s))]
  b <- ref$prob[match(keys, key(ref))]
  a[is.na(a)] <- 0
  b[is.na(b)] <- 0
  hour <- sub(" [0-9]+$", "", keys)
  tapply(abs(a - b), hour, sum) / 2
}
res <- do.call(rbind, lapply(c(1, 20, 64), function(nodes) {
  tv <- tv_per_hour(full(nodes))
  data.frame(nodes = nodes, animal_hours = length(tv), max_tv = max(tv), median_tv = median(tv))
}))
print(res, digits = 3)
write.csv(res, "validation/outputs/01_node_error.csv", row.names = FALSE)

# ---- Session ----
sessionInfo()
