# Title: Time of location_surface() for one animal of the shipped example
# Description: The numbers that ?location_surface states under Cost. On the
#   shipped example (stoney_track, 70 receivers, 4,818 water cells) at a
#   receiver-hour spread of 2.55 logits and the default 64 nodes, one hour
#   and all detected hours of stoney-01 are timed five times each, and the
#   medians are printed. Elapsed time depends on the machine and its load;
#   the load is printed with the result.
# Run from the package root: Rscript validation/03_location_surface_timing.R
# Output: the printed lines below, kept in
#   validation/outputs/03_location_surface_timing.txt.

# ---- Libraries ----
library(detfield)

# ---- Field and counts of the example ----
lake <- terra::rast(system.file("extdata", "stoney_depth.tif", package = "detfield"))
truth <- stoney_track$truth
curve <- de_logistic(d50 = truth$d50, d05 = truth$d05, depth = truth$depth)
field <- de_field(water_distances(lake, stoney_track$receivers),
                  de_model(curve, sigma = list(station_hour = 2.55)))
counts <- bin_detections(stoney_track$detections, stoney_track$receivers, stoney_track$tags)
one <- counts[counts$animal_id == "stoney-01", ]
hours <- unique(one$bin_start)
length(hours)  # detected hours of stoney-01
nrow(one)  # receiver-hour rows of stoney-01

# ---- Five timed runs ----
# A first call loads the code paths, so the timed calls are warm.
invisible(location_surface(one[one$bin_start == hours[1], ], field))
runs <- 5
times <- t(vapply(seq_len(runs), function(r) {
  c(one_hour = system.time(location_surface(one[one$bin_start == hours[1], ],
                                            field))[["elapsed"]],
    all_hours = system.time(location_surface(one, field))[["elapsed"]])
}, numeric(2)))
print(times)  # elapsed seconds, one row per run
cat(sprintf("stoney-01, median of %d runs: 1 hour %.2f s, %d hours %.2f s\n",
            runs, stats::median(times[, "one_hour"]), length(hours),
            stats::median(times[, "all_hours"])))
parallel::detectCores()  # cores of this machine
if (nzchar(Sys.which("uptime"))) system("uptime", intern = TRUE)  # load at the end

# ---- Session ----
sessionInfo()
