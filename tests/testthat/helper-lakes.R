# A 1,200 m by 600 m lake on a 30 m grid, with a peninsula down from the
# north shore in column 20, rows 1 to 14.
peninsula_lake <- function() {
  grid <- terra::rast(nrows = 20, ncols = 40, xmin = 0, xmax = 1200, ymin = 0, ymax = 600,
                      crs = "EPSG:32617", vals = 1)
  grid[terra::cellFromRowColCombine(grid, 1:14, 20)] <- NA
  grid
}

# Receivers at cell centres either side of the peninsula: (row 6, col 10)
# and (row 6, col 30).
peninsula_receivers <- function() {
  data.frame(station_id = c("west", "east"), x = c(285, 885), y = c(435, 435),
             deploy_z = c(4, 6), stringsAsFactors = FALSE)
}

pair <- function(d, id, cell) d[d[[attr(d, "key")]] == id & d$cell == cell, ]
