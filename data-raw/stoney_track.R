# ---- Build the Stoney Lake example ----
# Run from the package root: Rscript data-raw/stoney_track.R
# Copies the Stoney depth raster to inst/extdata/stoney_depth.tif, places a
# dense positioning array (70 receivers on a grid), walks five animals for 48
# hours, transmits and detects under the Stoney range model, and writes
# data/stoney_track.rda and the files inst/extdata/stoney_*.csv. The build
# itself is build_stoney_track() in data-raw/stoney_track_build.R. The range
# model is static: no variation between receivers or hours, so a surface
# computed with that curve and a spread of 0 uses the model the detections
# came from.

# ---- Libraries ----
library(terra)
pkgload::load_all(".", quiet = TRUE)
source("data-raw/stoney_track_build.R")

# ---- Settings ----
source_tif <- "../arraydesignR/inst/extdata/stoney_depth.tif"
tif_md5 <- "21ba8e62d4f1fe7050b155a8fa6c59b3"
start <- as.POSIXct("2024-06-01 00:00:00", tz = "UTC")

# ---- Depth raster ----
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
file.copy(source_tif, "inst/extdata/stoney_depth.tif", overwrite = TRUE)
md5 <- unname(tools::md5sum("inst/extdata/stoney_depth.tif"))
cat("md5:", md5, "\n")
stopifnot(md5 == tif_md5)
grid <- rast("inst/extdata/stoney_depth.tif")
grid
# Depth is stored as metres of zero or below; abs() gives depth below the surface
cat("water cells:", sum(!is.na(values(grid))), "of", ncell(grid), "\n")
range(abs(values(grid)), na.rm = TRUE)

# ---- Build the study ----
stoney_track <- build_stoney_track(grid, n_receivers = 70, n_animals = 5, hours = 48,
                                   start = start)
nrow(stoney_track$receivers)  # receivers
range(stoney_track$receivers$deploy_z)  # receiver depth, m
nrow(stoney_track$detections)  # detections, all animals
table(stoney_track$detections$animal_id)  # detections per animal
table(stoney_track$transmissions$animal_id)  # transmissions sent per animal
# Transmissions heard by at least one receiver, per animal
tapply(stoney_track$transmissions$n_detections > 0, stoney_track$transmissions$animal_id, sum)
table(stoney_track$track$animal_id)  # minute positions per animal
# Path length per animal, km
trk <- stoney_track$track
round(tapply(seq_len(nrow(trk)), trk$animal_id, function(i) {
  sum(sqrt(diff(trk$x[i])^2 + diff(trk$y[i])^2))
}) / 1000, 1)

# ---- Save the dataset ----
save(stoney_track, file = "data/stoney_track.rda", compress = "xz")
cat("data/stoney_track.rda:", file.size("data/stoney_track.rda"), "bytes\n")

# ---- Write the files ----
# The same study as files of the kind a study's own export gives: times as
# UTC text to the whole second, column names of CONTRACT.md.
csv_time <- function(x) {
  for (col in names(x)) if (inherits(x[[col]], "POSIXct")) {
    x[[col]] <- format(x[[col]], "%Y-%m-%d %H:%M:%S", tz = "UTC")
  }
  x
}
deploy_cols <- c("receiver_serial", "station_id", "deploy_datetime", "recover_datetime",
                 "deploy_lat", "deploy_lon", "deploy_z")
utils::write.csv(csv_time(stoney_track$receivers[, deploy_cols]),
                 "inst/extdata/stoney_deployments.csv", row.names = FALSE)
utils::write.csv(csv_time(stoney_track$detections), "inst/extdata/stoney_detections.csv",
                 row.names = FALSE)
utils::write.csv(stoney_track$tags, "inst/extdata/stoney_tags.csv", row.names = FALSE)
print(file.size(list.files("inst/extdata", pattern = "^stoney", full.names = TRUE)))
