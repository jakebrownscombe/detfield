# Code of the vignette "Hourly positions against the true track".
# Extracted from vignettes/positioning.Rmd by tools/purl_vignettes.R, a repository file
# the installed package does not carry. The vignette is the source; this file follows it.

# ---- Field ----
# Attach detfield, and ggplot2 for the figures
library(detfield)
library(ggplot2)
# Read the lake raster and the receivers of the example study
lake <- terra::rast(system.file("extdata", "stoney_depth.tif",
                                package = "detfield"))
receivers <- stoney_track$receivers
# Build the field from the curve the detections were drawn from, spread 0
truth <- stoney_track$truth
curve <- de_logistic(d50 = truth$d50, d05 = truth$d05, depth = truth$depth)
field <- de_field(water_distances(lake, receivers),
                  de_model(curve, sigma = list(station_hour = 0)))

# ---- Surfaces ----
# Count detections per animal, hour and listening receiver
counts <- bin_detections(stoney_track$detections, receivers, stoney_track$tags)
# Compute one surface per animal and detected hour
surface <- location_surface(counts, field)
surface

# ---- Estimates ----
# Take each stored cell's position and probability
cells <- data.frame(animal_id = surface$animal_id,
                    bin_start = surface$bin_start, prob = surface$prob,
                    round(terra::xyFromCell(lake, surface$cell), 1))
# Average the cell positions of each hour, weighted by probability
estimates <- aggregate(cbind(x = x * prob, y = y * prob, prob) ~
                         animal_id + bin_start, data = cells, FUN = sum)
estimates[c("x", "y")] <- estimates[c("x", "y")] / estimates$prob
# Count the receivers that heard the animal in each hour
heard <- counts[counts$y > 0, ]
receivers_heard <- aggregate(list(receivers_heard = heard$station_id),
                             by = heard[c("animal_id", "bin_start")],
                             FUN = length)
estimates <- merge(estimates, receivers_heard)
table(estimates$receivers_heard)

# ---- Truth ----
# Give each true minute position its hour, leaving out the final midnight
track <- stoney_track$track
track <- track[track$datetime < max(track$datetime), ]
track$bin_start <- as.POSIXct(trunc(track$datetime, "hours"))
# Average the minute positions of each hour
true_hours <- aggregate(cbind(x_true = x, y_true = y) ~
                          animal_id + bin_start, data = track, FUN = mean)
# Join each estimate to its truth and measure the error, in metres
hours <- merge(estimates, true_hours)
hours$error_m <- sqrt((hours$x - hours$x_true)^2 + (hours$y - hours$y_true)^2)
summary(hours$error_m)

# ---- Window ----
# Turn the study's share of single-receiver hours into a count out of 12
n_single <- round(12 * mean(hours$receivers_heard == 1))
n_single
# Test each run of 12 hours from each detected hour, by animal and time
hours <- hours[order(hours$animal_id, hours$bin_start), ]
fits_rule <- function(i) {
  run <- hours[hours$animal_id == hours$animal_id[i] &
                 hours$bin_start >= hours$bin_start[i] &
                 hours$bin_start < hours$bin_start[i] + 12 * 3600, ]
  nrow(run) == 12 && sum(run$receivers_heard == 1) == n_single &&
    run$receivers_heard[1] != 1
}
first <- Find(fits_rule, seq_len(nrow(hours)))
animal <- hours$animal_id[first]
window <- hours$bin_start[first] + c(0, 12 * 3600)
in_window <- function(time) time >= window[1] & time < window[2]
# Select the animal's track, transmissions and hourly positions in the window
window_track <- track[track$animal_id == animal & in_window(track$datetime), ]
transmissions <- stoney_track$transmissions
window_tx <- transmissions[transmissions$animal_id == animal &
                             in_window(transmissions$datetime), ]
window_hours <- hours[hours$animal_id == animal & in_window(hours$bin_start), ]
window_hours[c("bin_start", "receivers_heard", "error_m")]

# ---- Window receivers ----
# Set the map window to the track and the positions, with 600 m round them
view <- list(x = range(window_track$x, window_hours$x) + c(-600, 600),
             y = range(window_track$y, window_hours$y) + c(-600, 600))
in_view <- function(points, margin = 150) {
  points$x > view$x[1] + margin & points$x < view$x[2] - margin &
    points$y > view$y[1] + margin & points$y < view$y[2] - margin
}
# Sum the animal's detections per receiver over the 12 hours
window_counts <- counts[counts$animal_id == animal &
                          in_window(counts$bin_start), ]
window_rx <- aggregate(list(detections = window_counts$y),
                       by = window_counts["station_id"], FUN = sum)
window_rx <- merge(receivers[c("station_id", "x", "y")], window_rx)
window_rx <- window_rx[in_view(window_rx), ]
# Turn depth to metres below the surface (the raster stores it as negative)
depth_cells <- round(terra::as.data.frame(abs(lake), xy = TRUE), 1)
names(depth_cells)[3] <- "depth_m"
c(heard = sum(window_rx$detections > 0), silent = sum(window_rx$detections == 0))

# ---- Result positions ----
# Map depth, track, receivers, then each estimate joined to its truth
single <- window_hours$receivers_heard == 1
# A dark ground for the keys of white marks
dark_key <- theme(legend.key = element_rect(fill = "grey45", colour = NA))
ggplot(depth_cells, aes(x, y)) +
  geom_raster(aes(fill = depth_m)) +
  scale_fill_depth() +
  geom_path(data = window_track, aes(colour = "True track"), linewidth = 0.5) +
  layer_receivers(window_rx, count = "detections") +
  geom_segment(data = window_hours, aes(xend = x_true, yend = y_true, colour = "Error"),
               linewidth = 0.7) +
  geom_point(data = window_hours, aes(x_true, y_true, colour = "Mean true position"),
             shape = 21, fill = "#2ECC40", size = 2.2) +
  geom_point(data = window_hours[!single, ], aes(colour = "Posterior mean"), shape = 21,
             fill = "white", size = 2.6) +
  geom_point(data = window_hours[single, ], aes(colour = "Single-receiver hour"),
             shape = 24, fill = "white", size = 2.8) +
  # A constant mapped to colour gives each mark an entry in the key
  scale_colour_manual(name = NULL,
                      values = c("True track" = "#2ECC40", "Mean true position" = "black",
                                 "Posterior mean" = "black", "Single-receiver hour" = "black",
                                 "Error" = "white"),
                      breaks = c("True track", "Mean true position", "Posterior mean",
                                 "Single-receiver hour", "Error"),
                      guide = guide_legend(order = 3, theme = dark_key)) +
  coord_equal(xlim = view$x, ylim = view$y) +
  labs(x = "Easting (m, UTM 17N)", y = "Northing (m, UTM 17N)") +
  annotation_scale_bar(location = "br") +
  theme_positionr()

# ---- Single receiver means ----
# Join each single-receiver hour to the receiver that heard the animal
single_hours <- merge(hours[hours$receivers_heard == 1, ],
                      heard[c("animal_id", "bin_start", "station_id")])
single_hours <- merge(single_hours, receivers[c("station_id", "x", "y")],
                      by = "station_id", suffixes = c("", "_receiver"))
# Measure the distance from the posterior mean to that receiver, in metres
single_hours$to_receiver_m <- sqrt((single_hours$x - single_hours$x_receiver)^2 +
                                     (single_hours$y - single_hours$y_receiver)^2)
round(summary(single_hours$to_receiver_m))

# ---- Track map ----
# Map depth, track, missed then heard transmissions, then receivers
ggplot(depth_cells, aes(x, y)) +
  geom_raster(aes(fill = depth_m)) +
  scale_fill_depth() +
  geom_path(data = window_track, aes(colour = "True track"), linewidth = 0.5) +
  geom_point(data = window_tx[window_tx$n_detections == 0, ],
             aes(colour = "Transmission no receiver heard"), shape = 4, alpha = 0.6,
             size = 1.3) +
  geom_point(data = window_tx[window_tx$n_detections > 0, ],
             aes(colour = "Transmission heard"), shape = 21, fill = "#2ECC40", size = 1.6,
             stroke = 0.3) +
  # A constant mapped to colour gives each mark an entry in the key
  scale_colour_manual(name = NULL,
                      values = c("True track" = "grey45", "Transmission heard" = "black",
                                 "Transmission no receiver heard" = "red3"),
                      breaks = c("True track", "Transmission heard",
                                 "Transmission no receiver heard"),
                      guide = guide_legend(order = 3)) +
  layer_receivers(window_rx, count = "detections") +
  coord_equal(xlim = view$x, ylim = view$y) +
  labs(x = "Easting (m, UTM 17N)", y = "Northing (m, UTM 17N)") +
  annotation_scale_bar(location = "br") +
  theme_positionr()

# ---- Hour panels data ----
# Select the window's surface cells and scale each hour to its largest cell
hour_cells <- cells[cells$animal_id == animal & in_window(cells$bin_start), ]
hour_cells$share <- hour_cells$prob /
  ave(hour_cells$prob, hour_cells$bin_start, FUN = max)
hour_cells <- hour_cells[hour_cells$share >= 0.02, ]
# Label each hour by its start and the receivers that heard the animal
label_hour <- function(bin_start) {
  n_rx <- window_hours$receivers_heard[match(bin_start, window_hours$bin_start)]
  sprintf("%s, %s", trimws(format(bin_start, "%e June %H:00", tz = "UTC")),
          ifelse(n_rx == 1, "single receiver", paste(n_rx, "receivers")))
}
panels <- label_hour(sort(window_hours$bin_start))
hour_cells$panel <- factor(label_hour(hour_cells$bin_start), panels)
window_track$panel <- factor(label_hour(window_track$bin_start), panels)
window_hours$panel <- factor(label_hour(window_hours$bin_start), panels)
# Give every water cell a row in every panel, empty outside the surface,
# so each panel is a full grid
water <- depth_cells[c("x", "y")]
hour_grid <- merge(water, data.frame(panel = factor(panels, panels)), by = NULL)
hour_grid <- merge(hour_grid, hour_cells[c("panel", "x", "y", "share")],
                   all.x = TRUE)
levels(hour_grid$panel)

# ---- Hour panels receivers ----
# List every receiver once per hour, with that hour's count of the animal
window_counts$panel <- factor(label_hour(window_counts$bin_start), panels)
hour_rx <- data.frame(panel = window_counts$panel,
                      station_id = window_counts$station_id,
                      detections = window_counts$y)
hour_rx <- merge(hour_rx, receivers[c("station_id", "x", "y")])
hour_rx <- hour_rx[in_view(hour_rx), ]
nrow(hour_rx)

# ---- Hour panels ----
# Map water, each hour's surface, its track segment, truth, estimate, rings
ggplot(hour_grid, aes(x, y)) +
  geom_raster(aes(fill = share)) +
  scale_fill_probability(name = "Probability,\nshare of the\nhour's largest",
                         na_value = "grey88") +
  geom_path(data = window_track, aes(colour = "Track in the hour"), linewidth = 0.6) +
  geom_point(data = window_hours, aes(x_true, y_true, colour = "Mean true position"),
             shape = 21, fill = "#2ECC40", size = 1.8) +
  geom_point(data = window_hours, aes(colour = "Posterior mean"), shape = 3, size = 2.2,
             stroke = 1) +
  # A constant mapped to colour gives each mark an entry in the key
  scale_colour_manual(name = NULL,
                      values = c("Track in the hour" = "#2ECC40",
                                 "Mean true position" = "black",
                                 "Posterior mean" = "white"),
                      breaks = c("Track in the hour", "Mean true position",
                                 "Posterior mean"),
                      guide = guide_legend(order = 3, theme = dark_key)) +
  layer_receivers(hour_rx, count = "detections", size = 1.4, max_size = 3.6,
                  stroke = 0.5) +
  facet_wrap(~panel, ncol = 3) +
  coord_equal(xlim = view$x, ylim = view$y) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  labs(x = NULL, y = NULL) +
  annotation_scale_bar(location = "br", width_hint = 0.3, text_size = 2.2) +
  theme_positionr(base_size = 10)

# ---- Error table ----
# Group the hours by receivers heard, four or more together
hours$heard_by <- factor(pmin(hours$receivers_heard, 4), 1:4,
                         c("1", "2", "3", "4 or more"))
# Summarise hours and error per group, in metres
data.frame(hours = tapply(hours$error_m, hours$heard_by, length),
           median_error_m = round(tapply(hours$error_m, hours$heard_by, median)),
           q90_error_m = round(tapply(hours$error_m, hours$heard_by, quantile,
                                      probs = 0.9)))
# Summarise each animal: hours, share heard by one receiver, median error
data.frame(hours = tapply(hours$error_m, hours$animal_id, length),
           single_receiver = round(tapply(hours$receivers_heard == 1,
                                          hours$animal_id, mean), 2),
           median_error_m = round(tapply(hours$error_m, hours$animal_id, median)))

# ---- Error plot ----
# Plot error against receivers heard, one point per hour, with the medians
ggplot(hours, aes(heard_by, error_m)) +
  geom_point(aes(colour = animal_id), alpha = 0.7, size = 1.6,
             position = position_jitter(width = 0.25, height = 0, seed = 1)) +
  stat_summary(aes(linetype = "Median of the group"), fun = median, fun.min = median,
               fun.max = median, geom = "crossbar", width = 0.6, linewidth = 0.4) +
  scale_colour_viridis_d(name = "Animal", end = 0.9) +
  scale_linetype_manual(name = NULL, values = "solid") +
  labs(x = "Receivers that heard the animal in the hour",
       y = "Error of the hourly position (m)") +
  theme_positionr()

# ---- Coverage ----
# Find the 95% region of every hour
regions <- credible_region(surface, level = 0.95)
region_area <- attr(regions, "summary")
region_key <- paste(regions$animal_id, regions$bin_start, regions$cell)
# Mark each true minute of a detected hour as inside its hour's region or not
track$cell <- terra::cellFromXY(lake, as.matrix(track[c("x", "y")]))
minutes <- merge(track[c("animal_id", "bin_start", "cell")],
                 hours[c("animal_id", "bin_start", "heard_by")])
minutes$inside <- paste(minutes$animal_id, minutes$bin_start,
                        minutes$cell) %in% region_key
# Mark each hour's mean true position as inside its own region or not
hours$true_cell <- terra::cellFromXY(lake,
                                     as.matrix(hours[c("x_true", "y_true")]))
hours$mean_inside <- paste(hours$animal_id, hours$bin_start,
                           hours$true_cell) %in% region_key
# Summarise per hour, then by receivers heard
per_hour <- aggregate(inside ~ animal_id + bin_start, data = minutes, FUN = mean)
per_hour <- merge(merge(per_hour, hours), region_area)
coverage <- data.frame(
  hours = tapply(per_hour$inside, per_hour$heard_by, length),
  area_ha = tapply(per_hour$area_m2 / 1e4, per_hour$heard_by, median),
  minutes_inside = round(tapply(per_hour$inside, per_hour$heard_by, median), 2),
  mean_inside = round(tapply(per_hour$mean_inside, per_hour$heard_by, mean), 2)
)
coverage
c(minutes_inside = round(mean(minutes$inside), 3),
  hours_mean_inside = sum(hours$mean_inside), hours = nrow(hours))

# ---- Use tables ----
# Average the hourly surfaces per cell over all detected hours
use <- aggregate(list(share = surface$prob), by = list(cell = surface$cell),
                 FUN = sum)
use$share <- use$share / nrow(hours)
# Share of the true minutes of the detected hours in each cell
true_use <- aggregate(list(share = minutes$cell),
                      by = list(cell = minutes$cell), FUN = length)
true_use$share <- true_use$share / nrow(minutes)
# Measure the overlap of the two: 1 when identical, 0 when disjoint
both <- merge(use, true_use, by = "cell", all = TRUE,
              suffixes = c("_estimated", "_true"))
both[is.na(both)] <- 0
overlap <- sum(pmin(both$share_estimated, both$share_true))
round(overlap, 2)

# ---- Use near track ----
# Distance from every cell to the nearest cell the animals swam through
visited <- terra::rast(lake)
visited[unique(minutes$cell)] <- 1
to_track <- terra::values(terra::distance(visited), mat = FALSE)
# Share of the use map within 200 m, and in the cells swum through
near_track <- c(within_200_m = sum(use$share[to_track[use$cell] <= 200]),
                same_cell = sum(use$share[to_track[use$cell] == 0]))
round(near_track, 2)

# ---- Use cut ----
# Keep, per map, the smallest set of cells holding 95% of its total
top_95 <- function(map) {
  map <- map[order(-map$share), ]
  map[cumsum(map$share) - map$share < 0.95 * sum(map$share), ]
}
use_95 <- top_95(use)
true_use_95 <- top_95(true_use)
c(estimated_cells = nrow(use_95), of = nrow(use),
  true_cells = nrow(true_use_95), of = nrow(true_use))

# ---- Use maps ----
# Stack the two maps as panels, each cell with its coordinates
panel_names <- c("(a) Estimated: mean of the hourly surfaces",
                 "(b) True: minutes of the detected hours")
use_cells <- rbind(data.frame(use_95, panel = panel_names[1]),
                   data.frame(true_use_95, panel = panel_names[2]))
use_cells <- cbind(use_cells, round(terra::xyFromCell(lake, use_cells$cell), 1))
# Give every water cell a row in both panels, empty where it has no share
use_cells <- merge(merge(water, data.frame(panel = panel_names), by = NULL),
                   use_cells, all.x = TRUE)
use_cells$panel <- factor(use_cells$panel, panel_names)
# Map water in grey, then the shares on one magma scale for both panels
ggplot(use_cells, aes(x, y)) +
  geom_raster(aes(fill = share)) +
  scale_fill_viridis_c(name = "Share of\nanimal-hours\nper cell",
                       option = "magma", trans = "sqrt", na.value = "grey88") +
  layer_receivers(receivers, size = 1.2, stroke = 0.45) +
  facet_wrap(~panel, ncol = 1) +
  coord_equal() +
  labs(x = "Easting (m, UTM 17N)", y = "Northing (m, UTM 17N)") +
  annotation_scale_bar(location = "br") +
  theme_positionr()

# ---- Use rings ----
# Split each drawn cell's probability by whether its hour had one receiver
cell_mass <- data.frame(animal_id = surface$animal_id,
                        bin_start = surface$bin_start, cell = surface$cell,
                        prob = surface$prob)
cell_mass <- merge(cell_mass[cell_mass$cell %in% use_95$cell, ],
                   hours[c("animal_id", "bin_start", "receivers_heard")])
cell_mass$single_prob <- cell_mass$prob * (cell_mass$receivers_heard == 1)
by_cell <- aggregate(cbind(prob, single_prob) ~ cell, data = cell_mass,
                     FUN = sum)
# Share of the drawn probability from single-receiver hours
round(sum(by_cell$single_prob) / sum(by_cell$prob), 2)
# Cells where single-receiver hours give most of it, by animal, in km
ring_cells <- by_cell$cell[by_cell$single_prob > by_cell$prob / 2]
rings <- aggregate(single_prob ~ animal_id + cell, FUN = sum,
                   data = cell_mass[cell_mass$cell %in% ring_cells, ])
rings <- rings[order(-rings$single_prob), ]
rings <- rings[!duplicated(rings$cell), ]
rings <- cbind(rings, terra::xyFromCell(lake, rings$cell) / 1000)
span <- function(v) paste(round(range(v), 1), collapse = " to ")
data.frame(cells = as.vector(table(rings$animal_id)),
           east_km = tapply(rings$x, rings$animal_id, span),
           north_km = tapply(rings$y, rings$animal_id, span))

# ---- Presence draws ----
# Draw 50 cells per detected hour in proportion to its surface
presence <- sample_locations(surface, n = 50, seed = 1)
attr(presence, "estimand")
nrow(presence)

# ---- Draws sf ----
# Turn the presence draws to points in the raster's EPSG code
draw_points <- sf::st_as_sf(presence, coords = c("x", "y"),
                            crs = attr(presence, "epsg"), remove = FALSE)
draw_points[1:3, c("animal_id", "bin_start", "draw", "cell")]

# ---- Write files ----
# Write the draws to a CSV and to a GeoPackage layer
csv_file <- file.path(tempdir(), "stoney_draws.csv")
utils::write.csv(sf::st_drop_geometry(draw_points), csv_file,
                 row.names = FALSE)
gpkg_file <- file.path(tempdir(), "stoney_positions.gpkg")
sf::st_write(draw_points, gpkg_file, layer = "draws",
             delete_dsn = file.exists(gpkg_file), quiet = TRUE)
sf::st_layers(gpkg_file)
