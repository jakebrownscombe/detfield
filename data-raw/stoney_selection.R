# ---- Build the Stoney habitat-selection example ----
# Run from the package root: Rscript data-raw/stoney_selection.R
# Forty animals walk the correlated random walk of
# arraydesign::crw_model(step_mean = 30, step_sd = 15, turn_sd = 40,
# time_step_s = 60) on Stoney Lake for 48 h with a known preference for bottom depth (arraydesign's
# habitat_selection()): selection-01 to -20 at beta = -0.3 per m (shallow),
# selection-21 to -40 at beta = 0 (none). They are heard by the 70 receivers
# of stoney_track under its range model. These are the trial's animals at
# those two betas on that array (same build functions and seeds), so the
# article's live numbers equal rows of the trial, which is included beside them.
# Included beside the study, so the article starts from it and builds no field:
# the position draws of the 40 animals (20 per detected hour, seed = the
# animal's replicate), the coverage of the array per day (system_de(), one
# transmission and 24), the trial's per-animal rows for two blocks of 20
# animals per beta (block 2 is a replication with its own animals) and the
# trial's hourly depths on 70 receivers.
# Needs arraydesign and the outputs of the trial, a simulation study run
# outside the package and not included in the repository. Point
# DETFIELD_TRIAL_DIR at the folder that holds them. data-raw/ is not part of
# the package.
# Then run data-raw/stoney_selection_track_xy.R, which adds the minute
# positions (x, y) to `track`.

# ---- Libraries ----
pkgload::load_all(".", quiet = TRUE)
source("data-raw/stoney_selection_build.R")
set <- selection_settings
trial_dir <- Sys.getenv("DETFIELD_TRIAL_DIR")
if (!nzchar(trial_dir) || !dir.exists(trial_dir)) {
  stop("Set DETFIELD_TRIAL_DIR to the folder that holds the trial outputs.", call. = FALSE)
}

# ---- Lake, array, field ----
lake <- terra::rast("inst/extdata/stoney_depth.tif")
starts <- selection_starts(lake)
receivers <- selection_receivers(lake, 70)
stopifnot(isTRUE(all.equal(receivers[c("x", "y")], stoney_track$receivers[c("x", "y")],
                           check.attributes = FALSE)))
cat("receivers:", nrow(receivers), "(the array of stoney_track)\n")
field <- selection_field(lake, receivers)
end_time <- set$start + set$hours * 3600

# ---- The two groups ----
groups <- list(list(beta = -0.3, offset = 0L), list(beta = 0, offset = 20L))
parts <- lapply(groups, function(g) {
  tracks <- selection_tracks(lake, g$beta, starts)
  k <- g$offset + seq_len(set$n_animals)
  study <- selection_study(tracks, field, receivers, sprintf("selection-%02d", k),
                           sprintf("A69-9001-%d", 3000 + k), g$beta)
  # The same detections as the trial's study at this beta
  trial <- readRDS(file.path(trial_dir, "outputs", "trial", "studies",
                             sprintf("dense_beta%+.1f.rds", g$beta)))
  stopifnot(identical(study$detections[c("datetime", "receiver_serial")],
                      trial$detections[c("datetime", "receiver_serial")]))
  study$detections <- study$detections[study$detections$datetime < end_time, ]
  study
})
detections <- do.call(rbind, lapply(parts, `[[`, "detections"))
detections <- detections[order(detections$datetime, detections$receiver_serial), ]
rownames(detections) <- NULL
tags <- do.call(rbind, lapply(parts, `[[`, "tags"))
track <- do.call(rbind, lapply(parts, `[[`, "track"))
rownames(track) <- NULL
cat("detections:", nrow(detections), "; animals heard:", length(unique(detections$animal_id)),
    "of", nrow(tags), "\n")
print(table(beta = tags$beta[match(detections$animal_id, tags$animal_id)]))
cat("minute positions:", nrow(track), "(40 animals x 2,881)\n")

# ---- Position draws of the 40 animals ----
# The trial's draws: 20 per detected hour from each hour's surface, with the
# animal's replicate number as seed.
counts <- bin_detections(detections, receivers, tags)
surface <- location_surface(counts, field)
presence <- do.call(rbind, lapply(seq_len(nrow(tags)), function(k) {
  hours <- surface[surface$animal_id == tags$animal_id[k], ]
  sample_locations(hours, n = 20, seed = tags$replicate[k])
}))
rownames(presence) <- NULL
attr(presence, "estimand") <- attr(sample_locations(surface[surface$animal_id ==
                                                              "selection-01", ], n = 1,
                                                    seed = 1), "estimand")
attr(presence, "epsg") <- attr(sample_locations(surface[1, ], n = 1, seed = 1), "epsg")
cat("presence draws:", nrow(presence), "(",
    length(unique(paste(presence$animal_id, presence$bin_start))),
    "detected hours x 20 )\n")

# ---- Coverage of the array, per day ----
listening <- expand.grid(station_id = receivers$station_id, date = as.Date(set$start) + 0:2,
                         stringsAsFactors = FALSE)
coverage <- list(
  transmission = system_de(field, listening)[c("date", "cell", "p")],
  hour = system_de(field, listening,
                   transmissions = 3600 / mean(set$tag_delay))[c("date", "cell", "p")]
)
cat("coverage rows per table:", nrow(coverage$hour), "(3 days x", nrow(coverage$hour) / 3,
    "water cells)\n")

# ---- The trial rows, two blocks ----
# Detected hours, both arrays, four availabilities, both sources, every beta:
# one row per animal.
rows <- readRDS(file.path(trial_dir, "outputs", "trial_blocks", "rows.rds"))
trial <- rows[rows$hours == "detected", c("block", "beta", "array", "availability", "source",
                                          "replicate", "estimate", "se", "fits_warned",
                                          "n_hours", "single_share")]
trial$block <- as.integer(sub("block", "", trial$block))
trial$array <- ifelse(trial$array == "dense", 70L, 20L)
names(trial)[names(trial) == "array"] <- "receivers"
names(trial)[names(trial) == "replicate"] <- "animal"
trial$warned <- trial$fits_warned > 0
trial$fits_warned <- NULL
trial <- trial[order(trial$block, trial$receivers, trial$availability, trial$source, trial$beta,
                     trial$animal), ]
rownames(trial) <- NULL
cat("trial rows:", nrow(trial), "\n")
print(table(block = trial$block, receivers = trial$receivers))

# ---- Hourly depths, 70 receivers, two blocks ----
hours <- readRDS(file.path(trial_dir, "outputs", "trial_blocks", "hours.rds"))
hours <- hours[hours$array == "dense", c("block", "beta", "replicate", "bin_start", "receivers",
                                         "true_m", "drawn_m")]
hours$block <- as.integer(sub("block", "", hours$block))
hours$receivers <- as.integer(hours$receivers)
names(hours)[names(hours) == "replicate"] <- "animal"
hours <- hours[order(hours$block, hours$beta, hours$animal, hours$bin_start), ]
rownames(hours) <- NULL
cat("hourly depth rows:", nrow(hours), "\n")
# The example's hours carry the depth of the draws included above
ex <- hours[hours$block == 1 & hours$beta %in% c(-0.3, 0), ]
lake_depth <- abs(terra::values(lake, mat = FALSE))
ids <- tags$animal_id[match(paste(ex$beta, ex$animal), paste(tags$beta, tags$replicate))]
drawn <- tapply(lake_depth[presence$cell], paste(presence$animal_id,
                                                 as.numeric(presence$bin_start)), mean)
stopifnot(nrow(ex) == length(drawn),
          max(abs(ex$drawn_m - drawn[paste(ids, as.numeric(ex$bin_start))])) < 1e-12)

# ---- Save ----
stoney_selection <- list(
  receivers = stoney_track$receivers, detections = detections, tags = tags, track = track,
  presence = presence, coverage = coverage,
  truth = list(layer = "bottom depth, m (abs of stoney_depth.tif)", betas = c(-0.3, 0),
               hours = set$hours, start = set$start, step_mean = set$step_mean,
               step_sd = set$step_sd, turn_sd = set$turn_sd, time_step = set$time_step,
               tag_delay = set$tag_delay, transmissions_per_hour = 3600 / mean(set$tag_delay),
               draws = 20L, available = 20L, seed_presence = "the animal's replicate number",
               seed_available = "100 + the animal's replicate number", seeds = set$seeds,
               seeds_block2 = c(receivers = 1L, starts = 2020L, walk = 1001L,
                                transmissions = 1002L, detections = 1003L)),
  trial = trial, hours = hours
)
save(stoney_selection, file = "data/stoney_selection.rda", compress = "xz")
cat("data/stoney_selection.rda:", file.size("data/stoney_selection.rda"), "bytes\n")
