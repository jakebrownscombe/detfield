# Title: The grid path of water_distances() against the exact path
# Description: The numbers that NEWS.md states about
#   control = list(path = "grid") and about the 8-neighbour grid path on
#   open water. On the shipped example (stoney_depth.tif, 100 m cells, and
#   the 70 receivers of stoney_track, no max_distance), water_distances() is
#   run with the exact path (the default) and the grid path, and the two
#   tables are compared pair by pair.
#   The grid option keeps the straight line on clear pairs and takes the
#   8-neighbour path of terra::gridDist() where land lies on the line.
#   positionR 1.0.0 and 1.1.0 took the 8-neighbour grid path for every pair,
#   from the centre of the receiver's cell; from 1.2.0 the straight line is
#   used where no land crosses it. Part 3 measures the 8-neighbour path on
#   the pairs whose straight line stays in water.
# Run from the package root:
#   Rscript -e 'source("validation/05_grid_path.R", echo = TRUE)' \
#     > validation/outputs/05_grid_path.log 2>&1
# Output: validation/outputs/05_grid_path.log

# ---- Libraries ----
library(detfield)

# ---- Both paths on the example ----
lake <- terra::rast(system.file("extdata", "stoney_depth.tif", package = "detfield"))
rx <- stoney_track$receivers
exact <- suppressWarnings(water_distances(lake, rx))
grid <- suppressWarnings(water_distances(lake, rx, control = list(path = "grid")))
nrow(exact)  # receiver-cell pairs, exact path
nrow(grid)  # receiver-cell pairs, grid option
pairs <- merge(as.data.frame(exact), as.data.frame(grid), by = c("station_id", "cell"),
               suffixes = c("_exact", "_grid"))
nrow(pairs)  # pairs in both tables
identical(pairs$in_water_exact, pairs$in_water_grid)  # same line-of-sight test

# ---- Part 1: pairs whose straight line stays in water ----
clear <- pairs[pairs$in_water_exact, ]
nrow(clear)  # clear pairs
sum(clear$distance_m_grid == clear$straight_m_grid)  # grid option gives the straight line
sum(clear$distance_m_exact == clear$straight_m_exact)  # exact path gives the straight line

# ---- Part 2: pairs with land on the straight line ----
blocked <- pairs[!pairs$in_water_exact, ]
ratio <- blocked$distance_m_grid / blocked$distance_m_exact
length(ratio)  # blocked pairs
range(ratio)  # grid option over exact path
mean(ratio)
stats::median(ratio)
sum(ratio < 0.99)  # grid option shorter by more than 1%
sum(ratio < 0.90)  # grid option shorter by more than 10%
sum(ratio > 1.01)  # grid option longer by more than 1%
# Pairs the grid option reaches and the exact path does not: water joined
# to the receiver's water only where two land cells meet diagonally.
nrow(grid) - nrow(pairs)
cat(sprintf("Blocked pairs: grid option over exact %.3f to %.3f, mean %.3f, over %s pairs\n",
            min(ratio), max(ratio), mean(ratio), format(length(ratio), big.mark = ",")))

# ---- Part 3: the 8-neighbour path on clear pairs ----
# Distance along the 8-neighbour grid from the centre of each receiver's cell,
# against the straight line between the two cell centres.
info <- detfield:::grid_info(lake)
cells <- attr(exact, "cells")
packed <- attr(exact, "grid")
rx_used <- attr(exact, "receivers")
xy <- terra::xyFromCell(lake, cells)
over <- unlist(lapply(seq_len(nrow(rx_used)), function(i) {
  path <- detfield:::grid_path(packed, info, cells, rx_used$cell[i], 1000)
  centre <- terra::xyFromCell(lake, rx_used$cell[i])
  straight <- sqrt((xy[, 1] - centre[1])^2 + (xy[, 2] - centre[2])^2)
  mine <- clear[clear$station_id == rx_used$station_id[i], ]
  at <- match(mine$cell, cells)
  keep <- straight[at] > 0
  path[at][keep] / straight[at][keep]
}))
length(over)  # clear pairs off the receiver's own cell
range(over)  # 8-neighbour path over the straight line
cat(sprintf("8-neighbour path on clear pairs: up to %.1f%% longer, mean %.1f%%, over %s pairs\n",
            100 * (max(over) - 1), 100 * (mean(over) - 1),
            format(length(over), big.mark = ",")))
# The bound for any 8-neighbour path: at 22.5 degrees off an axis it is
# sqrt(4 - 2 * sqrt(2)) times the straight line.
sqrt(4 - 2 * sqrt(2))

# ---- Session ----
sessionInfo()
