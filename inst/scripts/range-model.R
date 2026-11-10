# Code of the vignette "Range models: bring your own, or start from two distances".
# Extracted from vignettes/range-model.Rmd by tools/purl_vignettes.R, a repository file
# the installed package does not carry. The vignette is the source; this file follows it.

# ---- Stand in ----
# Attach detfield, and ggplot2 for the figures
library(detfield)
library(ggplot2)
# Build the curve from d50 and d05 at 1 m and 35 m, with no spread
stand_in <- de_model(de_logistic(d50 = c(250, 750), d05 = c(500, 1500),
                                 depth = c(1, 35)),
                     sigma = list(station_hour = 0))
stand_in

# ---- Plot stand in ----
# Draw the curve at three receiver depths
plot(stand_in, depth = c(1, 10, 25))

# ---- Depth model ----
# Fit heard out of sent on distance, receiver depth and their interaction
depth_fit <- glm(cbind(heard, sent - heard) ~ distance_m * depth_m,
                 family = binomial, data = stoney_range_test)
depth_model <- de_model(depth_fit, sigma = list(station_hour = 0))
# Compare it with the stand-in at 250 m and 500 m, at 1 m and 25 m depth
check <- expand.grid(distance_m = c(250, 500), depth_m = c(1, 25))
check$fitted <- plogis(depth_model$predict(check$distance_m,
                                           depth_m = check$depth_m))
check$stand_in <- plogis(stand_in$predict(check$distance_m,
                                          depth_m = check$depth_m))
round(check, 3)

# ---- Prepare ----
# Number the days from the first day of the test
range_test <- stoney_range_covariates
first_day <- as.Date("2024-05-22")
range_test$day <- as.numeric(as.Date(range_test$datetime) - first_day)
# Make factors for the receiver and receiver-hour random effects
range_test$receiver_hour <- factor(paste(range_test$receiver_serial,
                                         format(range_test$datetime)))
range_test$receiver_serial <- factor(range_test$receiver_serial)
range_test$tag_model <- factor(range_test$tag_model)
c(rows = nrow(range_test), days = length(unique(range_test$day)),
  receiver_hours = nlevels(range_test$receiver_hour))

# ---- Fit time ----
# Fit a binomial GAM with a smooth over day and two random effects
time_fit <- mgcv::bam(
  cbind(heard, sent - heard) ~ distance_m * depth_m + tag_model +
    s(day, k = 5) + s(receiver_serial, bs = "re") +
    s(receiver_hour, bs = "re"),
  family = binomial, data = range_test, discrete = TRUE
)
summary(time_fit)$s.table

# ---- Check day ----
# Predict the curve's logit at 300 m and 1 m depth for a V13-1x tag, each day
days <- data.frame(day = 0:5, distance_m = 300, depth_m = 1,
                   tag_model = "V13-1x",
                   receiver_serial = levels(range_test$receiver_serial)[1],
                   receiver_hour = levels(range_test$receiver_hour)[1])
fitted_logit <- predict(time_fit, days,
                        exclude = c("s(receiver_serial)", "s(receiver_hour)"))
# Set the day effect the counts were simulated with, from the day's wind
wind <- tapply(range_test$wind_ms, range_test$day, unique)
simulated_logit <- -0.25 * (wind - 4)
# Centre both on their mean and compare
round(data.frame(day = 0:5, wind_ms = wind,
                 fitted = fitted_logit - mean(fitted_logit),
                 simulated = simulated_logit - mean(simulated_logit)), 2)

# ---- Spread ----
# Take the standard deviation of each random effect: for a binomial model
# its variance is 1 over the effect's smoothing parameter
sds <- 1 / sqrt(time_fit$sp[c("s(receiver_serial)", "s(receiver_hour)")])
round(sds, 2)
spread <- sqrt(sum(sds^2))
spread

# ---- Time model ----
# Turn a date into the day number the smooth was fitted on
day_of <- function(date) as.numeric(as.Date(date) - first_day)
# Add a receiver and a receiver-hour level, which exclude = leaves out
add_levels <- function(pairs) {
  pairs$receiver_serial <- levels(range_test$receiver_serial)[1]
  pairs$receiver_hour <- levels(range_test$receiver_hour)[1]
  pairs
}
time_model <- de_model(time_fit, sigma = list(station_hour = spread),
                       modulators = list(day = day_of),
                       covariates = add_levels,
                       exclude = c("s(receiver_serial)", "s(receiver_hour)"))
time_model

# ---- Distances ----
# Read the lake and measure distances through water from the receivers
lake <- terra::rast(system.file("extdata", "stoney_depth.tif",
                                package = "detfield"))
receivers <- stoney_track$receivers
distances <- water_distances(lake, receivers)

# ---- Time fields ----
# Build one field per day, for a V13-1x tag
field_days <- as.Date(c("2024-05-22", "2024-05-24"))
day_fields <- lapply(field_days, function(date) {
  suppressMessages(de_field(distances, time_model, tag_model = "V13-1x",
                            date = date))
})
# Share of water cells at least one receiver hears at 0.5 or more, per day
any_receiver <- lapply(day_fields, function(field) {
  1 - apply(1 - field$p, 1, prod)
})
round(vapply(any_receiver, function(p) mean(p >= 0.5), 0), 3)

# ---- Plot time fields ----
# Put both days' probabilities in one table, one panel per day
panels <- paste0(c("(a) ", "(b) "), format(field_days, "%e %B"), ", wind ",
                 wind[c("0", "2")], " m/s")
day_cells <- do.call(rbind, Map(function(field, p, panel) {
  data.frame(terra::xyFromCell(lake, field$cells), p = p, panel = panel)
}, day_fields, any_receiver, panels))
day_cells$panel <- factor(day_cells$panel, panels)
# Map both, with the receivers
ggplot(day_cells, aes(x, y)) +
  geom_raster(aes(fill = p)) +
  scale_fill_probability(name = "Detection\nprobability") +
  layer_receivers(receivers, size = 0.9, stroke = 0.35) +
  facet_wrap(~panel, nrow = 1) +
  coord_equal() +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  labs(x = NULL, y = NULL) +
  annotation_scale_bar(location = "br", text_size = 2.4) +
  theme_positionr()

# ---- Listening ----
# Build the field of the stand-in curve
field <- de_field(distances, stand_in)
# Mark all 70 receivers listening on 1 June, and all but S034 on 2 June
stations <- receivers$station_id
listening <- data.frame(
  station_id = c(stations, setdiff(stations, "S034")),
  date = as.Date(rep(c("2024-06-01", "2024-06-02"),
                     c(length(stations), length(stations) - 1)))
)
table(listening$date)

# ---- Coverage ----
# Compute coverage by at least one, and at least two, listening receivers
cover_1 <- system_de(field, listening, k = 1)
cover_2 <- system_de(field, listening, k = 2)
# Share of water cells at 0.5 or more, per day
daily_shares <- data.frame(
  date = sort(unique(listening$date)),
  heard_by_one = tapply(cover_1$p >= 0.5, cover_1$date, mean),
  heard_by_two = tapply(cover_2$p >= 0.5, cover_2$date, mean)
)
daily_shares

# ---- Coverage map data ----
# Give each day's coverage its coordinates and a panel label
day_titles <- c("(a) 1 June, all receivers", "(b) 2 June, S034 not listening")
cover_1 <- cbind(cover_1, terra::xyFromCell(lake, cover_1$cell))
cover_1$panel <- factor(day_titles[match(cover_1$date,
                                         sort(unique(cover_1$date)))],
                        day_titles)
# Set a 4 km window round S034, and mark it on the second day
s034 <- receivers[receivers$station_id == "S034", ]
view <- list(x = s034$x + c(-2200, 2200), y = s034$y + c(-2200, 2000))
lost <- data.frame(s034[c("x", "y")], panel = factor(day_titles[2],
                                                     day_titles))
nrow(cover_1)

# ---- Coverage map ----
# Map the coverage, the listening receivers and a cross at S034 on day two
ggplot(cover_1, aes(x, y)) +
  geom_raster(aes(fill = p)) +
  scale_fill_probability(name = "Detection\nprobability") +
  layer_receivers(receivers[receivers$station_id != "S034", ], size = 2) +
  layer_receivers(transform(s034, panel = factor(day_titles[1], day_titles)),
                  size = 2) +
  geom_point(data = lost, shape = 4, colour = "#1A1A1A", size = 3,
             stroke = 2.4) +
  geom_point(data = lost, shape = 4, colour = "#FDE725", size = 3,
             stroke = 1.2) +
  facet_wrap(~panel, nrow = 1) +
  coord_equal(xlim = view$x, ylim = view$y) +
  labs(x = "Easting (m, UTM 17N)", y = "Northing (m, UTM 17N)") +
  annotation_scale_bar(location = "br") +
  theme_positionr()

# ---- Silent hour ----
# Count the detections of the example study per animal, hour and receiver
counts <- bin_detections(stoney_track$detections, receivers,
                         stoney_track$tags)
# Find the busiest receiver-hour, and the counts of that animal-hour
busiest <- counts[which.max(counts$y), ]
as_logged <- counts[counts$animal_id == busiest$animal_id &
                      counts$bin_start == busiest$bin_start, ]
as_logged[as_logged$y > 0, c("animal_id", "bin_start", "station_id", "y", "n")]

# ---- Silent cases ----
# Set the busiest receiver's count to 0, still marked listening
silent_by_mistake <- as_logged
silent_by_mistake$y[silent_by_mistake$station_id == busiest$station_id] <- 0L
# Leave the busiest receiver out of the hour
left_out <- as_logged[as_logged$station_id != busiest$station_id, ]
hour_versions <- list(as_logged = as_logged,
                      silent_by_mistake = silent_by_mistake,
                      left_out = left_out)

# ---- Silent receiver ----
# Find the cells of the animal's true positions in the hour
track <- stoney_track$track
in_hour <- track$animal_id == busiest$animal_id &
  track$datetime >= busiest$bin_start &
  track$datetime < busiest$bin_start + 3600
true_cells <- terra::cellFromXY(lake, as.matrix(track[in_hour, c("x", "y")]))
# Compute each version's surface and 95% region, and score it
silent_scores <- t(vapply(hour_versions, function(hour_counts) {
  region <- credible_region(location_surface(hour_counts, field))
  c(region_ha = attr(region, "summary")$area_m2 / 1e4,
    minutes_inside = sum(true_cells %in% region$cell),
    minutes = length(true_cells))
}, numeric(3)))
silent_scores
