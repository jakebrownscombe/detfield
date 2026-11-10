# ---- Build the Stoney Lake range tests ----
# Run from the package root: Rscript data-raw/stoney_range_test.R
# Writes the files of a simulated range test to inst/extdata
# (stoney_range_test_receivers.csv, stoney_range_test_tags.csv,
# stoney_range_test_detections.csv) and two counted range tests to data/.
# stoney_range_test is what range_test_counts() makes of those files: it
# follows the Stoney range model (d50 250 m and d05 500 m at 1 m receiver
# depth, 750 m and 1,500 m at 35 m) at receiver depths of 1 m and 25 m,
# which bracket the 20 receivers of the array (1.0 to 24.8 m), with no
# variation between receivers or hours.
# stoney_range_covariates adds two tag models, daily wind and a spread
# between receivers and hours, for fitting a range model with covariates.

pkgload::load_all(".", quiet = TRUE)

# ---- Settings ----
set.seed(20261102)
utc <- function(x) as.POSIXct(x, tz = "UTC")
stoney_curve <- de_logistic(d50 = c(250, 750), d05 = c(500, 1500), depth = c(1, 35))
tag_delay <- c(120, 180)
test_depth <- c(1, 25)

# ---- stoney_range_test: the layout ----
# Two lines, one in 1 m of water and one in 25 m. On each line three
# receivers sit at 0, 450 and 900 m and five range-test tags at 150, 350,
# 650, 1,100 and 1,600 m from the line's west end, so the 15 pairs of a line
# span 100 to 1,600 m. Positions are laid out on an azimuthal equidistant
# projection centred on each line's west end, the stations S006 (1.0 m) and
# S016 (24.8 m) of the 20-receiver grid array of arraydesign's examples,
# and written as WGS 84. The shallow line ran from 19 to 22 May and the deep
# line from 22 to 25 May, so no receiver of one line was in the water with a
# tag of the other.
load("data/stoney_track.rda")
ends <- stoney_track$receivers[match(c("S006", "S016"), stoney_track$receivers$station_id), ]
lon_lat <- function(east_m, end) {
  aeqd <- sprintf("+proj=aeqd +lat_0=%f +lon_0=%f +datum=WGS84 +units=m", end$deploy_lat,
                  end$deploy_lon)
  pts <- sf::st_as_sf(data.frame(x = east_m, y = 0), coords = c("x", "y"), crs = aeqd)
  round(sf::st_coordinates(sf::st_transform(pts, 4326)), 7)
}
rx_at <- c(0, 450, 900)
tag_at <- c(150, 350, 650, 1100, 1600)
shift <- function(k) (k - 1) * 3 * 86400
rt_receivers <- do.call(rbind, lapply(1:2, function(k) {
  ll <- lon_lat(rx_at, ends[k, ])
  data.frame(receiver_serial = as.character(480100 + 10 * k + 1:3),
             station_id = paste0("RT", k, 1:3),
             deploy_datetime = utc("2024-05-19 23:40:00") + shift(k),
             recover_datetime = utc("2024-05-22 00:20:00") + shift(k),
             deploy_lat = ll[, "Y"], deploy_lon = ll[, "X"], deploy_z = test_depth[k])
}))
rt_tags <- do.call(rbind, lapply(1:2, function(k) {
  ll <- lon_lat(tag_at, ends[k, ])
  data.frame(transmitter = sprintf("A69-9001-%d", 9000 + 10 * k + 1:5), tag_model = "V13-1x",
             deploy_datetime = utc("2024-05-19 23:50:00") + shift(k),
             recover_datetime = utc("2024-05-22 00:10:00") + shift(k),
             deploy_lat = ll[, "Y"], deploy_lon = ll[, "X"],
             min_delay = tag_delay[1], max_delay = tag_delay[2])
}))
rt_receivers
rt_tags

# ---- stoney_range_test: the transmissions and detections ----
# Each tag sends from its deployment with a delay drawn uniformly from 120 to
# 180 s.
# Each receiver on the tag's line hears each transmission with the Stoney curve's
# probability at their distance and the line's depth, independently.
sends <- do.call(rbind, lapply(seq_len(nrow(rt_tags)), function(t) {
  first <- as.numeric(rt_tags$deploy_datetime[t])
  last <- as.numeric(rt_tags$recover_datetime[t])
  at <- first + runif(1, 0, tag_delay[2]) + cumsum(c(0, runif(1500, tag_delay[1], tag_delay[2])))
  data.frame(transmitter = rt_tags$transmitter[t], line = (t - 1) %/% 5 + 1,
             east = tag_at[(t - 1) %% 5 + 1], time = at[at < last])
}))
nrow(sends)  # transmissions, all ten tags
heard <- do.call(rbind, lapply(seq_len(nrow(rt_receivers)), function(r) {
  k <- (r - 1) %/% 3 + 1
  mine <- sends[sends$line == k, ]
  d <- abs(rx_at[(r - 1) %% 3 + 1] - mine$east)
  p <- plogis(stoney_curve$predict(d, depth_m = test_depth[k]))
  keep <- runif(nrow(mine)) < p
  data.frame(datetime = as.POSIXct(floor(mine$time[keep]), origin = "1970-01-01", tz = "UTC"),
             receiver_serial = rt_receivers$receiver_serial[r],
             transmitter = mine$transmitter[keep])
}))
rt_detections <- heard[order(heard$datetime, heard$receiver_serial), ]
rownames(rt_detections) <- NULL
nrow(rt_detections)  # detections, all six receivers

# ---- stoney_range_test: the files counted ----
stoney_range_test <- range_test_counts(rt_detections, rt_receivers, rt_tags)
nrow(stoney_range_test)  # receiver, tag and hour rows
print(aggregate(cbind(heard, sent) ~ depth_m + distance_m, stoney_range_test, sum))

csv_time <- function(x) {
  for (col in names(x)) if (inherits(x[[col]], "POSIXct")) {
    x[[col]] <- format(x[[col]], "%Y-%m-%d %H:%M:%S", tz = "UTC")
  }
  x
}
utils::write.csv(csv_time(rt_receivers), "inst/extdata/stoney_range_test_receivers.csv",
                 row.names = FALSE)
utils::write.csv(csv_time(rt_tags), "inst/extdata/stoney_range_test_tags.csv",
                 row.names = FALSE)
utils::write.csv(csv_time(rt_detections), "inst/extdata/stoney_range_test_detections.csv",
                 row.names = FALSE)

# ---- stoney_range_covariates ----
set.seed(20261103)
# Six receivers, three at 1 m and three at 25 m, and eight range-test tags of
# two models, each model at four positions. Six days of hours; one wind speed
# per day. The V13-1x curve at a wind of 4 m/s and no receiver or hour effect
# is the Stoney curve at each depth. A V9-2x tag is 0.9 logits lower, each
# m/s of wind above 4 lowers the logit by 0.25, each receiver adds a normal
# effect with SD 0.5, centred within its depth, and each receiver-hour one
# with SD 0.9.
rx <- data.frame(receiver_serial = as.character(480201:480206),
                 station_id = paste0("RD", 1:6), depth_m = rep(test_depth, each = 3),
                 rx_pos = rep(c(0, 400, 800), 2))
tags <- data.frame(transmitter = sprintf("A69-9001-%d", 9101:9108),
                   tag_model = rep(c("V9-2x", "V13-1x"), each = 4),
                   tag_pos = rep(c(150, 400, 800, 1400), 2))
pairs <- merge(rx, tags, by = NULL)
pairs$distance_m <- abs(pairs$rx_pos - pairs$tag_pos)
hours <- seq(utc("2024-05-22 00:00:00"), by = 3600, length.out = 6 * 24)
wind <- data.frame(date = as.Date("2024-05-22") + 0:5, wind_ms = c(2, 5, 8, 3, 6, 4))
grid <- merge(pairs, data.frame(datetime = hours), by = NULL)
grid$wind_ms <- wind$wind_ms[match(as.Date(grid$datetime), wind$date)]
# Only three receivers sit at each depth, so the receiver effects are centred
# within each depth; otherwise their mean at one depth moves that depth's curve.
receiver_effect <- rnorm(6, 0, 0.5)
receiver_effect <- setNames(receiver_effect - ave(receiver_effect, rx$depth_m), rx$receiver_serial)
receiver_effect
key <- paste(grid$receiver_serial, as.numeric(grid$datetime))
hour_effect <- rnorm(length(unique(key)), 0, 0.9)[match(key, unique(key))]
eta <- stoney_curve$predict(grid$distance_m, depth_m = grid$depth_m) -
  0.9 * (grid$tag_model == "V9-2x") - 0.25 * (grid$wind_ms - 4) +
  receiver_effect[grid$receiver_serial] + hour_effect
# 24 transmissions an hour: 3,600 s over the mean delay of 150 s.
grid$sent <- 24L
grid$heard <- rbinom(nrow(grid), grid$sent, plogis(eta))
stoney_range_covariates <- grid[order(grid$datetime, grid$receiver_serial, grid$transmitter),
                                c("datetime", "receiver_serial", "station_id", "depth_m",
                                  "transmitter", "tag_model", "distance_m", "wind_ms",
                                  "sent", "heard")]
rownames(stoney_range_covariates) <- NULL
nrow(stoney_range_covariates)  # receiver, tag and hour rows

# ---- Save ----
save(stoney_range_test, file = "data/stoney_range_test.rda", compress = "xz")
save(stoney_range_covariates, file = "data/stoney_range_covariates.rda", compress = "xz")
cat("data/stoney_range_test.rda:", file.size("data/stoney_range_test.rda"), "bytes\n")
cat("data/stoney_range_covariates.rda:", file.size("data/stoney_range_covariates.rda"),
    "bytes\n")
