# ---- Add the true positions to the stoney_selection track ----
# Run from the package root, after data-raw/stoney_selection.R:
#   Rscript data-raw/stoney_selection_track_xy.R
# stoney_selection.R kept only the cell of each minute position. A track drawn
# through cell centres steps along the 100 m grid and, over 48 h, retraces the
# same cell edges, so on a map it reads as a lattice rather than a path. The
# simulator's own positions are regenerated here from the same settings and
# seeds, checked cell by cell against the shipped track, and added as `x` and
# `y` (metres, UTM 17N). Each is kept to 10 m from its cell centre, cut toward
# the centre, so every position stays in its own cell. Plain rounding would
# move positions near a cell edge into the next cell. 10 m keeps the package
# tarball under the 2 MB gate: at 0.1 m the data file was 874,512 bytes and the
# tarball 2,206,323; at 10 m the file is about 617,500. The cells, which the
# analyses read, are exact. Nothing else in the object changes.
# Needs arraydesign; data-raw/ is not part of the package.

# ---- Libraries ----
pkgload::load_all(".", quiet = TRUE)
source("data-raw/stoney_selection_build.R")

# ---- Regenerate the tracks of the two groups ----
lake <- terra::rast("inst/extdata/stoney_depth.tif")
starts <- selection_starts(lake)
groups <- list(list(beta = -0.3, offset = 0L), list(beta = 0, offset = 20L))
tracks <- do.call(rbind, lapply(groups, function(g) {
  tr <- selection_tracks(lake, g$beta, starts)
  # Offsets from the cell centre, cut to 10 m toward the centre
  centre <- terra::xyFromCell(lake, tr$cell)
  data.frame(animal_id = sprintf("selection-%02d", g$offset + tr$path_id),
             datetime = tr$time, cell = as.integer(tr$cell),
             x = centre[, 1] + trunc((tr$x - centre[, 1]) / 10) * 10,
             y = centre[, 2] + trunc((tr$y - centre[, 2]) / 10) * 10)
}))
nrow(tracks)  # minute positions, 40 animals x 2,881

# ---- Check against the shipped track ----
study <- stoney_selection
stopifnot(identical(tracks$animal_id, study$track$animal_id),
          identical(tracks$datetime, study$track$datetime),
          identical(tracks$cell, study$track$cell))
cat("cells match the shipped track at all", nrow(tracks), "minute positions\n")
# Each position falls in its own cell
in_cell <- terra::cellFromXY(lake, as.matrix(tracks[c("x", "y")])) == tracks$cell
cat("positions inside their cell:", sum(in_cell), "of", length(in_cell), "\n")
stopifnot(all(in_cell))

# ---- Save ----
study$track <- tracks
stoney_selection <- study
usethis::use_data(stoney_selection, overwrite = TRUE, compress = "xz")
cat("data/stoney_selection.rda:", file.size("data/stoney_selection.rda"), "bytes\n")
