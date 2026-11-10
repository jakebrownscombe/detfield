# ---- Shortest paths through water ----

random_lake <- function(seed, n = 30, land = 0.15, res = 30) {
  set.seed(seed)
  grid <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n * res, ymin = 0, ymax = n * res,
                      crs = "EPSG:32617", vals = 1)
  grid[sample(terra::ncell(grid), round(land * terra::ncell(grid)))] <- NA
  grid
}

test_that("the path around one island is the closed-form shortest path", {
  # 30 m cells; one land cell at x 150 to 180, y 150 to 180. The line from
  # (15, 165) to (315, 165) runs through it, so the path bends at the corners
  # (150, 180) and (180, 180), or the two below: 2 * sqrt(135^2 + 15^2) + 30.
  grid <- terra::rast(nrows = 11, ncols = 11, xmin = 0, xmax = 330, ymin = 0, ymax = 330,
                      crs = "EPSG:32617", vals = 1)
  grid[terra::cellFromXY(grid, cbind(165, 165))] <- NA
  d <- water_distances(grid, data.frame(station_id = "w", x = 15, y = 165))
  across <- d[d$cell == terra::cellFromXY(grid, cbind(315, 165)), ]
  expect_false(across$in_water)
  expect_equal(across$distance_m, 2 * sqrt(135^2 + 15^2) + 30, tolerance = 1e-9)
  # Two cells up, the line passes above the island and is used.
  up <- d[d$cell == terra::cellFromXY(grid, cbind(315, 225)), ]
  expect_true(up$in_water)
  expect_equal(up$distance_m, sqrt(300^2 + 60^2))
})

test_that("distances between cell centres are symmetric", {
  grid <- random_lake(3)
  info <- grid_info(grid)
  at <- sample(which(info$water), 25)
  xy <- xy_of_cell(info, at)
  d <- water_distances(grid, data.frame(station_id = sprintf("S%02d", 1:25), x = xy[, 1],
                                        y = xy[, 2]))
  m <- matrix(NA_real_, 25, 25)
  key <- match(d$station_id, sprintf("S%02d", 1:25))
  k <- match(d$cell, at)
  m[cbind(key[!is.na(k)], k[!is.na(k)])] <- d$distance_m[!is.na(k)]
  expect_equal(m, t(m), tolerance = 1e-9)
})

test_that("the exact path agrees with an independent visibility graph", {
  # The reference shares no code with detfield: corners found here, line of
  # sight by GEOS on the water polygon, shortest paths by Floyd-Warshall. The
  # land rule closes the point where two land cells meet diagonally, so a
  # GEOS line through such a point counts as blocked.
  grid <- random_lake(5, n = 25)
  info <- grid_info(grid)
  polys <- terra::as.polygons(terra::ifel(is.na(grid), NA, 1))
  water <- sf::st_set_crs(sf::st_union(sf::st_as_sf(polys)), NA)
  land <- matrix(is.na(terra::values(grid, mat = FALSE)), info$nrow, byrow = TRUE)
  pad <- matrix(TRUE, info$nrow + 2, info$ncol + 2)
  pad[2:(info$nrow + 1), 2:(info$ncol + 1)] <- land
  g <- expand.grid(i = 0:info$nrow, j = 0:info$ncol)
  q <- cbind(pad[cbind(g$i + 1, g$j + 1)], pad[cbind(g$i + 1, g$j + 2)],
             pad[cbind(g$i + 2, g$j + 1)], pad[cbind(g$i + 2, g$j + 2)])
  to_xy <- function(k) cbind(info$xmin + g$j[k] * 30, info$ymax - g$i[k] * 30)
  corners <- to_xy(which(rowSums(q) == 1))
  pinch <- to_xy(which(rowSums(q) == 2 & ((q[, 1] & q[, 4]) | (q[, 2] & q[, 3]))))
  pinch_sf <- sf::st_sfc(lapply(seq_len(nrow(pinch)), function(i) sf::st_point(pinch[i, ])))
  sees <- function(a, b) {
    s <- sf::st_sfc(lapply(seq_len(nrow(a)), function(i) sf::st_linestring(rbind(a[i, ], b[i, ]))))
    ok <- as.vector(sf::st_covers(water, s, sparse = FALSE))
    if (length(pinch_sf)) ok <- ok & lengths(sf::st_is_within_distance(s, pinch_sf, 1e-3)) == 0
    ok
  }
  n <- nrow(corners)
  pr <- which(upper.tri(diag(n)), arr.ind = TRUE)
  dmat <- matrix(Inf, n, n)
  diag(dmat) <- 0
  vis <- sees(corners[pr[, 1], , drop = FALSE], corners[pr[, 2], , drop = FALSE])
  len <- sqrt(rowSums((corners[pr[, 1], , drop = FALSE] - corners[pr[, 2], , drop = FALSE])^2))
  dmat[pr[vis, , drop = FALSE]] <- len[vis]
  dmat[pr[vis, 2:1, drop = FALSE]] <- len[vis]
  for (k in seq_len(n)) dmat <- pmin(dmat, outer(dmat[, k], dmat[k, ], "+"))
  cells <- which(info$water)
  cxy <- xy_of_cell(info, cells)
  set.seed(9)
  rx <- cxy[sample(nrow(cxy), 8), ]
  d <- water_distances(grid, data.frame(station_id = sprintf("R%d", 1:8), x = rx[, 1], y = rx[, 2]))
  cc <- expand.grid(c = seq_len(nrow(cxy)), v = seq_len(n))
  cell_corner <- matrix(Inf, nrow(cxy), n)
  ok <- sees(cxy[cc$c, , drop = FALSE], corners[cc$v, , drop = FALSE])
  cell_corner[cbind(cc$c, cc$v)[ok, , drop = FALSE]] <-
    sqrt(rowSums((cxy[cc$c, , drop = FALSE] - corners[cc$v, , drop = FALSE])^2))[ok]
  ref <- unlist(lapply(1:8, function(r) {
    first <- rep(Inf, n)
    s <- sees(matrix(rx[r, ], n, 2, byrow = TRUE), corners)
    first[s] <- sqrt(colSums((t(corners) - rx[r, ])^2))[s]
    dv <- apply(first + dmat, 2, min)
    path <- apply(sweep(cell_corner, 2, dv, "+"), 1, min)
    direct <- sees(matrix(rx[r, ], nrow(cxy), 2, byrow = TRUE), cxy)
    ifelse(direct, sqrt(colSums((t(cxy) - rx[r, ])^2)), path)
  }))
  ours <- rep(NA_real_, length(ref))
  pos <- (match(d$station_id, sprintf("R%d", 1:8)) - 1) * nrow(cxy) + match(d$cell, cells)
  ours[pos] <- d$distance_m
  ref[!is.finite(ref)] <- NA
  expect_gt(sum(!is.na(ref)), 4000)
  expect_identical(is.na(ours), is.na(ref))
  expect_lt(max(abs(ours - ref), na.rm = TRUE), 1e-6)
})

test_that("moving receivers by millimetres moves distances by millimetres", {
  grid <- random_lake(8, n = 40)
  info <- grid_info(grid)
  set.seed(4)
  at <- xy_of_cell(info, sample(which(info$water), 15))
  # Receivers off their cell centres, as a station file places them
  rx <- data.frame(station_id = sprintf("S%02d", 1:15), x = at[, 1] + stats::runif(15, -14, 14),
                   y = at[, 2] + stats::runif(15, -14, 14))
  base <- water_distances(grid, rx)
  for (rep in 1:3) {
    moved <- rx
    moved$x <- moved$x + stats::runif(15, -0.009, 0.009)
    moved$y <- moved$y + stats::runif(15, -0.009, 0.009)
    d <- water_distances(grid, moved)
    m <- merge(base, d, by = c("station_id", "cell"))
    expect_equal(nrow(m), nrow(base))
    expect_lt(max(abs(m$distance_m.x - m$distance_m.y)), 0.0128)
  }
})

test_that("a receiver on the point where two land cells meet is not shut in", {
  # 30 m cells; land at x 30 to 60, y 60 to 90 and x 60 to 90, y 30 to 60,
  # meeting at (60, 60). A receiver 1.4 cm from that point, in the water cell
  # north-east of it, must not be shut in by the point it sits beside.
  grid <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 120, ymin = 0, ymax = 120,
                      crs = "EPSG:32617", vals = 1)
  grid[terra::cellFromRowColCombine(grid, 2, 2)] <- NA
  grid[terra::cellFromRowColCombine(grid, 3, 3)] <- NA
  info <- grid_info(grid)
  xy <- xy_of_cell(info, which(info$water))
  seen <- line_clear(info, 60.01, 60.01, xy[, 1], xy[, 2])
  expect_equal(sum(seen), 8)
  expect_true(seen[xy[, 1] == 75 & xy[, 2] == 75])
  # It sees through the gap it sits in, and not through land.
  expect_true(seen[xy[, 1] == 15 & xy[, 2] == 15])
  expect_false(seen[xy[, 1] == 15 & xy[, 2] == 75])
  # A line from farther away through the same point stays blocked.
  expect_false(line_clear(info, 75, 75, 15, 15))
  d <- water_distances(grid, data.frame(station_id = "gap", x = 60.01, y = 60.01))
  expect_equal(nrow(d), 14)
  expect_equal(attr(d, "unreachable"), 0)
})

test_that("water joined only where two land cells meet is unreachable", {
  grid <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 120, ymin = 0, ymax = 120,
                      crs = "EPSG:32617", vals = 1)
  # A pond in the bottom-right cell, closed off by land on its north and west
  # but for the corner it shares with the cell diagonally north-west of it.
  grid[terra::cellFromRowColCombine(grid, 3, 4)] <- NA
  grid[terra::cellFromRowColCombine(grid, 4, 3)] <- NA
  d <- water_distances(grid, data.frame(station_id = "a", x = 15, y = 105))
  pond <- terra::cellFromRowCol(grid, 4, 4)
  expect_false(pond %in% d$cell)
  expect_equal(attr(d, "unreachable"), 1)
  g <- water_distances(grid, data.frame(station_id = "a", x = 15, y = 105),
                       control = list(path = "grid"))
  expect_true(pond %in% g$cell)
})

test_that("several cores give the same exact paths as one", {
  skip_on_os("windows")
  grid <- random_lake(12, n = 30)
  info <- grid_info(grid)
  set.seed(2)
  at <- xy_of_cell(info, sample(which(info$water), 6))
  rx <- data.frame(station_id = sprintf("S%d", 1:6), x = at[, 1], y = at[, 2])
  one <- water_distances(grid, rx)
  two <- water_distances(grid, rx, control = list(cores = 2))
  expect_identical(as.data.frame(unclass(two)), as.data.frame(unclass(one)))
  near <- water_distances(grid, rx, max_distance = 300)
  m <- merge(one, near, by = c("station_id", "cell"))
  expect_equal(nrow(near), sum(one$distance_m <= 300))
  expect_equal(m$distance_m.x, m$distance_m.y)
})
