# Code of the vignette "Fish positioning with detection range information".
# Extracted from vignettes/detfield.Rmd by tools/purl_vignettes.R, a repository file
# the installed package does not carry. The vignette is the source; this file follows it.

# ---- Read range test ----
# Attach detfield, and ggplot2 for the figures
library(detfield)
library(ggplot2)
# Find a file included with the package
in_file <- function(name) system.file("extdata", name, package = "detfield")
serial_text <- c(receiver_serial = "character")
# Read the three range-test files and count their rows
rt_detections <- read.csv(in_file("stoney_range_test_detections.csv"),
                          colClasses = serial_text)
rt_receivers <- read.csv(in_file("stoney_range_test_receivers.csv"),
                         colClasses = serial_text)
rt_tags <- read.csv(in_file("stoney_range_test_tags.csv"))
c(detections = nrow(rt_detections), receivers = nrow(rt_receivers),
  tags = nrow(rt_tags))

# ---- Range test times ----
# Convert detection, deployment and recovery times to UTC date-times
rt_detections$datetime <- as.POSIXct(rt_detections$datetime, tz = "UTC")
for (col in c("deploy_datetime", "recover_datetime")) {
  rt_receivers[[col]] <- as.POSIXct(rt_receivers[[col]], tz = "UTC")
  rt_tags[[col]] <- as.POSIXct(rt_tags[[col]], tz = "UTC")
}
range(rt_detections$datetime)

# ---- Count range test ----
# Count transmissions sent and heard per receiver, tag and hour
range_test <- range_test_counts(rt_detections, rt_receivers, rt_tags)
head(range_test[, c("datetime", "depth_m", "distance_m", "sent", "heard")])
# Count the rows per receiver depth (m) and distance (m)
table(range_test$depth_m, range_test$distance_m)

# ---- Fit range model ----
# Fit a binomial GLM of heard out of sent on distance and receiver depth
fit <- glm(cbind(heard, sent - heard) ~ distance_m * depth_m,
           family = binomial, data = range_test)
# Compute the Pearson dispersion: near 1 means no spread beyond the binomial
dispersion <- sum(residuals(fit, type = "pearson")^2) / df.residual(fit)
dispersion

# ---- De model ----
# Wrap the fitted model with a spread of 0
model <- de_model(fit, sigma = list(station_hour = 0))
model

# ---- Plot range model ----
# Predict the fitted curve every 10 m at the two tested depths
curves <- expand.grid(distance_m = seq(0, 1600, by = 10), depth_m = c(1, 25))
curves$p <- plogis(model$predict(curves$distance_m, depth_m = curves$depth_m))
# Pool the share heard at each depth and distance over hours
pooled <- aggregate(cbind(heard, sent) ~ depth_m + distance_m,
                    data = range_test, FUN = sum)
# Plot the curves over the pooled shares
ggplot(curves, aes(distance_m, p, colour = factor(depth_m))) +
  geom_line(linewidth = 1) +
  geom_point(data = pooled, aes(y = heard / sent), shape = 21,
             fill = "white", size = 2.4, stroke = 1) +
  scale_colour_viridis_d(name = "Receiver depth (m)", end = 0.7) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Distance (m)", y = "Detection probability per transmission") +
  theme_positionr()

# ---- Read lake ----
# Read the depth raster of the lake
lake <- terra::rast(in_file("stoney_depth.tif"))
# Read the receiver deployments and convert their times to UTC
deployments <- read.csv(in_file("stoney_deployments.csv"),
                        colClasses = serial_text)
for (col in c("deploy_datetime", "recover_datetime")) {
  deployments[[col]] <- as.POSIXct(deployments[[col]], tz = "UTC")
}
nrow(deployments)

# ---- Water distances ----
# Measure the distance through water from each receiver to each water cell
distances <- water_distances(lake, deployments)
distances

# ---- Field ----
# Compute the detection probability per receiver and water cell
field <- de_field(distances, model)
field

# ---- Plot field ----
# Map the field: the probability that any receiver detects a transmission
plot(field)

# ---- Read detections ----
# Read the detections and convert their times to UTC
detections <- read.csv(in_file("stoney_detections.csv"),
                       colClasses = serial_text)
detections$datetime <- as.POSIXct(detections$datetime, tz = "UTC")
# Read the tags and their delays
tags <- read.csv(in_file("stoney_tags.csv"))
tags

# ---- Bin detections ----
# Count detections per animal, hour and listening receiver
counts <- bin_detections(detections, deployments, tags)
head(counts[counts$y > 0, c("animal_id", "bin_start", "station_id", "y", "n")])
# Count the detected hours per animal
table(unique(counts[c("animal_id", "bin_start")])$animal_id)

# ---- Location surface ----
# Compute one surface per animal and detected hour
surface <- location_surface(counts, field)
surface

# ---- One hour ----
# Select the counts, at the receivers that heard it, from 07:00 UTC
hour_start <- as.POSIXct("2024-06-01 07:00:00", tz = "UTC")
hour_counts <- counts[counts$animal_id == "stoney-01" &
                        counts$bin_start == hour_start, ]
hour_counts[hour_counts$y > 0, c("station_id", "y", "n")]
# Select the true path of the hour, one position a minute
true_path <- stoney_track$track
true_path <- true_path[true_path$animal_id == "stoney-01" &
                         true_path$datetime >= hour_start &
                         true_path$datetime <= hour_start + 3600, ]
nrow(true_path)

# ---- Plot surface ----
# Map the surface, with the true path drawn on a white edge
plot(surface, animal_id = "stoney-01", bin_start = hour_start) +
  geom_path(data = true_path, aes(x, y), colour = "white", linewidth = 1.6,
            inherit.aes = FALSE) +
  geom_path(data = true_path, aes(x, y), colour = "grey45", linewidth = 0.8,
            inherit.aes = FALSE)

# ---- Credible region ----
# Find the 95% region of every hour, and its area
regions <- credible_region(surface, level = 0.95)
region_areas <- attr(regions, "summary")
# Compute the median area of the 95% region per animal, in ha
aggregate(list(median_ha = region_areas$area_m2 / 1e4),
          by = region_areas["animal_id"], FUN = function(ha) round(median(ha)))

# ---- Sample locations ----
# Draw 20 cells per detected hour in proportion to probability
draws <- sample_locations(surface, n = 20, seed = 1)
head(draws, 3)
attr(draws, "estimand")

# ---- Plot draws ----
# Turn depth to metres below the surface (the raster stores it as negative)
depth_cells <- terra::as.data.frame(abs(lake), xy = TRUE)
names(depth_cells)[3] <- "depth_m"
# Warm colours, which read on the blue depth ramp, one per animal; each
# track is a darker shade of its animal's colour, so it can be followed
# through its own draws
animals <- sort(unique(draws$animal_id))
draw_colours <- setNames(c("#FF8C00", "#FF3B3B", "#FF66CC", "#FFE135",
                           "#B3FF00"), animals)
track_colours <- grDevices::rgb(t(grDevices::col2rgb(draw_colours)) * 0.55,
                                maxColorValue = 255)
names(track_colours) <- animals
# Draw each animal's track thin, in its own dark shade
track_layers <- lapply(animals, function(a) {
  geom_path(data = stoney_track$track[stoney_track$track$animal_id == a, ],
            colour = track_colours[[a]], linewidth = 0.25)
})
# Map the depth, the draws, the tracks over them, then the receivers
ggplot(depth_cells, aes(x, y)) +
  geom_raster(aes(fill = depth_m)) +
  scale_fill_depth() +
  geom_point(data = draws, aes(colour = animal_id), alpha = 0.6, size = 1.3,
             stroke = 0) +
  scale_colour_manual(name = "Animal", values = draw_colours) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.6))) +
  track_layers +
  layer_receivers(attr(distances, "receivers"), size = 1.2,
                  stroke = 0.45) +
  coord_equal() +
  labs(x = "Easting (m, UTM 17N)", y = "Northing (m, UTM 17N)") +
  annotation_scale_bar(location = "br") +
  theme_positionr()
