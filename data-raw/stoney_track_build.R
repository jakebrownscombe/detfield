# ---- Build function for the Stoney Lake example ----
# Sourced by data-raw/stoney_track.R. Uses the
# arraydesign simulator, which detfield does not import; data-raw/ is not
# part of the package.
#
# The dense positioning example: the receivers are arraydesign's
# place_receivers() with n = 70, method = "grid" and seed = 1, on cells at
# least 1 m deep; the range model is the Stoney curve with no variation
# between receivers or hours; the animals walk a correlated random walk from
# random water cells; tags send every 120 to
# 180 s. Seeds are fixed so the object is reproducible.

build_stoney_track <- function(grid, n_receivers = 70, n_animals = 5, hours = 48,
                               start = as.POSIXct("2024-06-01 00:00:00", tz = "UTC"),
                               seeds = c(receivers = 1L, tracks = 1L, transmissions = 2L,
                                         detections = 3L)) {
  # The Stoney range model: d50 and d05 at 1 m and at 35 m receiver depth
  range_depth <- c(1, 35)
  range_d50 <- c(250, 750)
  range_d05 <- c(500, 1500)
  # Movement: 30 m a minute, turns of SD 40 degrees. Tags: 120 to 180 s delay
  step_mean <- 30
  step_sd <- 15
  turn_sd <- 40
  time_step <- 60
  tag_delay <- c(120, 180)
  # Receivers go in 12 h before the first position and come out 12 h after the last
  deploy_at <- start - 12 * 3600
  recover_at <- start + (hours + 12) * 3600

  # ---- Receivers ----
  deployable <- terra::ifel(abs(grid) >= 1, grid, NA)
  layout <- arraydesign::place_receivers(deployable, n = n_receivers, method = "grid",
                                         seed = seeds[["receivers"]])
  receivers <- data.frame(
    receiver_serial = layout$receiver_serial,
    station_id = layout$station_id,
    deploy_datetime = deploy_at,
    recover_datetime = recover_at,
    deploy_lat = round(layout$deploy_lat, 6),
    deploy_lon = round(layout$deploy_lon, 6),
    deploy_z = round(layout$deploy_z, 1),
    x = layout$x, y = layout$y,
    stringsAsFactors = FALSE
  )
  attr(receivers, "epsg") <- 32617L

  # ---- Detection field of the range model ----
  model <- de_model(de_logistic(range_d50, range_d05, depth = range_depth),
                    sigma = list(station_hour = 0))
  field <- suppressMessages(de_field(water_distances(grid, receivers), model))

  # ---- Tracks ----
  walk <- arraydesign::crw_model(step_mean = step_mean, step_sd = step_sd,
                                 turn_sd = turn_sd, time_step_s = time_step)
  tracks <- arraydesign::simulate_tracks(
    grid, walk, n_animals = n_animals, n_steps = hours * 3600 / time_step,
    control = arraydesign::track_control(start_time = start), seed = seeds[["tracks"]]
  )
  ids <- sprintf("stoney-%02d", seq_len(n_animals))
  transmitters <- sprintf("A69-9001-%d", 1000 + seq_len(n_animals))

  # ---- Transmissions and detections ----
  tx <- arraydesign::simulate_transmissions(tracks, tag_delay = tag_delay,
                                            seed = seeds[["transmissions"]])
  sim <- arraydesign::simulate_detections(tx, field, overdispersion = "none",
                                          seed = seeds[["detections"]])
  animal_of <- function(path_id) ids[path_id]
  detections <- sim$detections
  # The simulator names animals by path number; give them the example's names
  path_of <- as.integer(detections$animal_id)
  detections$animal_id <- animal_of(path_of)
  detections$transmitter <- transmitters[path_of]
  detections <- detections[order(detections$datetime, detections$receiver_serial),
                           c("datetime", "receiver_serial", "transmitter", "animal_id")]
  rownames(detections) <- NULL
  transmissions <- sim$transmissions
  transmissions$animal_id <- animal_of(transmissions$path_id)
  transmissions$transmitter <- transmitters[transmissions$path_id]
  transmissions <- transmissions[, c("animal_id", "transmitter", "datetime", "x", "y", "cell",
                                     "n_detections")]
  rownames(transmissions) <- NULL

  # ---- Known tracks ----
  # Positions at full precision: rounding moves a position at the shore onto land
  track <- data.frame(animal_id = animal_of(tracks$path_id), datetime = tracks$time,
                      x = tracks$x, y = tracks$y)

  tags <- data.frame(transmitter = transmitters, animal_id = ids,
                     min_delay = tag_delay[1], max_delay = tag_delay[2])

  list(
    receivers = receivers,
    detections = detections,
    tags = tags,
    track = track,
    transmissions = transmissions,
    truth = list(
      d50 = range_d50, d05 = range_d05, depth = range_depth, station_hour = 0,
      step_mean = step_mean, step_sd = step_sd, turn_sd = turn_sd, time_step = time_step,
      tag_delay = tag_delay, seed = seeds
    )
  )
}
