# ---- Build functions for the habitat-selection study ----
# Sourced by data-raw/stoney_selection.R, which ships stoney_selection. The
# trial whose results ship in stoney_selection$trial used these functions
# and seeds. Uses the arraydesign simulator, which
# detfield does not import; data-raw/ is not part of the package.
#
# The Stoney example with a known preference for bottom depth: animals walk
# the correlated random walk of
# arraydesign::crw_model(step_mean = 30, step_sd = 15, turn_sd = 40,
# time_step_s = 60)
# with arraydesign's habitat_selection() on abs(lake), beta per m, for 48 h
# from 20 fixed water cells; tags send every 120 to 180 s; detections are
# drawn from the Stoney range model (no spread between receivers or hours)
# on the 70-receiver array of stoney_track or the 20-receiver grid array of
# arraydesign's examples. Seeds are fixed
# so the objects are reproducible; the replicate is the animal.

selection_settings <- list(
  betas = c(-0.3, -0.1, 0, 0.1, 0.3), n_animals = 20L, hours = 48L,
  start = as.POSIXct("2024-06-01 00:00:00", tz = "UTC"),
  step_mean = 30, step_sd = 15, turn_sd = 40, time_step = 60, tag_delay = c(120, 180),
  d50 = c(250, 750), d05 = c(500, 1500), range_depth = c(1, 35),
  seeds = c(receivers = 1L, starts = 20L, walk = 1L, transmissions = 2L, detections = 3L)
)

# ---- Starts ----
# The same 20 water cells at every beta and on both arrays, so animal i is
# paired across the design.
selection_starts <- function(lake, set = selection_settings) {
  wg <- arraydesign:::water_grid(lake)
  withr::with_seed(set$seeds[["starts"]],
                   arraydesign:::cell_xy(wg, sample(which(wg$wet), set$n_animals)))
}

# ---- Arrays ----
# n = 70 is stoney_track's array; n = 20 is the grid array of arraydesign's
# examples.
selection_receivers <- function(lake, n, set = selection_settings) {
  deployable <- terra::ifel(abs(lake) >= 1, lake, NA)
  lay <- arraydesign::place_receivers(deployable, n = n, method = "grid",
                                      seed = set$seeds[["receivers"]])
  start <- set$start
  data.frame(receiver_serial = lay$receiver_serial, station_id = lay$station_id,
             deploy_datetime = start - 12 * 3600,
             recover_datetime = start + (set$hours + 12) * 3600,
             deploy_lat = round(lay$deploy_lat, 6), deploy_lon = round(lay$deploy_lon, 6),
             deploy_z = round(lay$deploy_z, 1), x = lay$x, y = lay$y,
             stringsAsFactors = FALSE)
}

selection_field <- function(lake, receivers, set = selection_settings) {
  model <- detfield::de_model(detfield::de_logistic(set$d50, set$d05, depth = set$range_depth),
                              sigma = list(station_hour = 0))
  suppressMessages(detfield::de_field(detfield::water_distances(lake, receivers), model))
}

# ---- Tracks ----
# Minute positions of 20 animals at one beta, with the minute's cell.
selection_tracks <- function(lake, beta, starts, set = selection_settings) {
  depth <- abs(lake)
  names(depth) <- "depth_m"
  walk <- arraydesign::crw_model(step_mean = set$step_mean, step_sd = set$step_sd,
                                 turn_sd = set$turn_sd, time_step_s = set$time_step)
  model <- arraydesign::habitat_selection(walk, depth, beta = beta)
  arraydesign::simulate_tracks(
    lake, model, n_animals = set$n_animals, n_steps = set$hours * 3600 / set$time_step,
    control = arraydesign::track_control(start = starts, start_time = set$start),
    seed = set$seeds[["walk"]]
  )
}

# ---- Detections ----
# Transmissions are drawn once per beta, so both arrays hear the same
# transmissions; `ids` and `transmitters` name the 20 animals.
selection_detections <- function(tracks, field, ids, transmitters, set = selection_settings) {
  tx <- arraydesign::simulate_transmissions(tracks, tag_delay = set$tag_delay,
                                            seed = set$seeds[["transmissions"]])
  sim <- arraydesign::simulate_detections(tx, field, overdispersion = "none",
                                          seed = set$seeds[["detections"]])
  det <- sim$detections
  path <- as.integer(det$animal_id)
  det$animal_id <- ids[path]
  det$transmitter <- transmitters[path]
  det <- det[order(det$datetime, det$receiver_serial),
             c("datetime", "receiver_serial", "transmitter", "animal_id")]
  rownames(det) <- NULL
  det
}

# ---- One study ----
# Tracks at one beta heard by one array: the tables location_surface() and
# the truth need, with animals named `ids`.
selection_study <- function(tracks, field, receivers, ids, transmitters, beta,
                            set = selection_settings) {
  track <- data.frame(animal_id = ids[tracks$path_id], datetime = tracks$time,
                      cell = as.integer(tracks$cell))
  tags <- data.frame(transmitter = transmitters, animal_id = ids, beta = beta,
                     replicate = seq_along(ids), min_delay = set$tag_delay[1],
                     max_delay = set$tag_delay[2], stringsAsFactors = FALSE)
  list(receivers = receivers,
       detections = selection_detections(tracks, field, ids, transmitters, set),
       tags = tags, track = track)
}
