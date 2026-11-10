# ---- The land test near grid lines, corners and closed gaps ----
# The first four tests put lines and receivers on, or within millimetres of,
# grid lines, corners and closed gaps, where a land test most often fails.

lake_with <- function(nr, nc, land_rc, res = 100) {
  g <- terra::rast(nrows = nr, ncols = nc, xmin = 0, xmax = nc * res, ymin = 0,
                   ymax = nr * res, crs = "EPSG:32617", vals = 1)
  g[terra::cellFromRowCol(g, land_rc[, 1], land_rc[, 2])] <- NA
  g
}

test_that("a line a millimetre off a grid line is blocked by land on both sides of it", {
  # 100 m cells, 10 x 4. Row 5 is land in columns 2 and 3, which meet along x = 200.
  info <- grid_info(lake_with(10, 4, cbind(5, 2:3)))
  for (off in c(0, 1e-3, -1e-3, 0.05, 0.5)) {
    expect_false(line_clear(info, 200 + off, 950, 200, 50), label = paste("offset", off))
    expect_false(line_clear(info, 200, 50, 200 + off, 950), label = paste("reverse, offset", off))
  }
})

test_that("line of sight is symmetric for receivers near grid lines", {
  set.seed(9)
  wet <- matrix(stats::runif(400) > 0.25, 20)
  info <- grid_info(lake_with(20, 20, which(!wet, arr.ind = TRUE)))
  p <- cbind(round(stats::runif(300, 0, 2000) / 100) * 100 + stats::runif(300, -0.1, 0.1),
             stats::runif(300, 0, 2000))
  q <- cbind(stats::runif(300, 0, 2000),
             round(stats::runif(300, 0, 2000) / 100) * 100 + stats::runif(300, -0.1, 0.1))
  expect_identical(line_clear(info, p[, 1], p[, 2], q[, 1], q[, 2]),
                   line_clear(info, q[, 1], q[, 2], p[, 1], p[, 2]))
})

test_that("a receiver just outside the tolerance of a closed gap sees its own cell", {
  # Land at rows 2 and 3, columns 2 and 3 diagonally: a closed gap at (200, 200).
  info <- grid_info(lake_with(4, 4, rbind(c(2, 2), c(3, 3))))
  for (a in c(10, 30, 45, 60, 80) * pi / 180) {
    for (r in c(0.05, 0.11, 0.13, 0.15, 0.5)) {
      x0 <- 200 + r * cos(a)
      y0 <- 200 + r * sin(a)
      expect_true(line_clear(info, x0, y0, 250, 250),
                  label = sprintf("angle %.0f, %.2f m", a * 180 / pi, r))
    }
  }
})

test_that("distances move by at most the move for receivers on grid lines and corners", {
  set.seed(5)
  wet <- matrix(stats::runif(400) > 0.2, 20)
  g <- lake_with(20, 20, which(!wet, arr.ind = TRUE))
  info <- grid_info(g)
  xy <- xy_of_cell(info, which(info$water))
  k <- sample(nrow(xy), 30)
  rx <- data.frame(station_id = sprintf("S%02d", 1:30),
                   x = xy[k, 1] + c(rep(50, 10), rep(49.999, 10), rep(49.95, 10)),
                   y = xy[k, 2] + c(rep(0, 10), rep(-49.999, 10), rep(10, 10)))
  rx$x <- pmin(rx$x, 1999.9)
  base <- suppressWarnings(water_distances(g, rx))
  moved <- rx
  a <- stats::runif(30, 0, 2 * pi)
  moved$x <- moved$x + 0.002 * cos(a)
  moved$y <- moved$y + 0.002 * sin(a)
  d <- suppressWarnings(water_distances(g, moved))
  ok <- attr(base, "receivers")$snap_m == 0 & attr(d, "receivers")$snap_m == 0
  m <- merge(base, d, by = c("station_id", "cell"))
  m <- m[m$station_id %in% rx$station_id[ok], ]
  expect_lt(max(abs(m$distance_m.x - m$distance_m.y)), 0.003)
})

test_that("receivers on and beside grid lines, corners, shores and closed gaps move by millimetres", { # nolint: line_length_linter.
  # Every placed point of placed_points(), moved 1 mm in a random direction.
  # A receiver that snaps to a water cell centre (it sat on land) is left out:
  # snapping moves it 50 to 70.7 m by design.
  g <- feature_lake()
  pts <- placed_points()
  base <- suppressWarnings(water_distances(g, pts[, c("station_id", "x", "y")]))
  set.seed(31)
  a <- stats::runif(nrow(pts), 0, 2 * pi)
  moved <- pts
  moved$x <- moved$x + 1e-3 * cos(a)
  moved$y <- moved$y + 1e-3 * sin(a)
  d <- suppressWarnings(water_distances(g, moved[, c("station_id", "x", "y")]))
  ok <- attr(base, "receivers")$snap_m == 0 & attr(d, "receivers")$snap_m == 0
  expect_gt(sum(ok), 200)
  m <- merge(base, d, by = c("station_id", "cell"), all = TRUE)
  m <- m[m$station_id %in% pts$station_id[ok], ]
  expect_false(anyNA(m$distance_m.x) || anyNA(m$distance_m.y))
  expect_lt(max(abs(m$distance_m.x - m$distance_m.y)), 1.001e-3 + 1e-9)
})

test_that("placed receivers get the distances of an independent reference", {
  # At least 200 receivers on and within 1 mm, 1 cm, 11 cm and 0.5 m of grid
  # lines, cell corners, shore corners and a closed gap, against
  # ref_distances() (helper-reference-path.R), at the position each was used at.
  g <- feature_lake()
  pts <- placed_points()
  expect_gte(nrow(pts), 200)
  d <- suppressWarnings(water_distances(g, pts[, c("station_id", "x", "y")]))
  rx <- attr(d, "receivers")
  cells <- attr(d, "cells")
  info <- grid_info(g)
  tgt <- xy_of_cell(info, cells)
  ref <- ref_distances(ref_build(g), cbind(rx$x, rx$y), tgt)
  ours <- matrix(Inf, nrow(rx), length(cells))
  ours[cbind(match(d$station_id, rx$station_id), match(d$cell, cells))] <- d$distance_m
  expect_identical(is.finite(ours), is.finite(ref))
  expect_lt(max(abs(ours - ref)[is.finite(ref)]), 1e-6)
})

test_that("de_model() on a discrete bam predicts each pair on its own", {
  skip_if_not_installed("mgcv")
  set.seed(1)
  n <- 5000
  dat <- data.frame(distance_m = stats::runif(n, 0, 1500), depth_m = stats::runif(n, 1, 30))
  dat$y <- stats::rbinom(n, 1, stats::plogis(2 - 0.004 * dat$distance_m))
  fit <- mgcv::bam(y ~ s(distance_m) + s(depth_m), family = binomial, data = dat, discrete = TRUE)
  m <- de_model(fit)
  batch <- c(300, stats::runif(4999, 0, 2000))
  expect_equal(m$predict(300, depth_m = 5), m$predict(batch, depth_m = 5)[1], tolerance = 1e-10)
  # de_model() passes discrete = FALSE to every bam; a fit that was not
  # discrete must predict as predict() does.
  plain <- mgcv::bam(y ~ s(distance_m) + s(depth_m), family = binomial, data = dat)
  nd <- data.frame(distance_m = batch[1:50], depth_m = 5)
  expect_equal(de_model(plain)$predict(batch[1:50], depth_m = 5),
               as.vector(stats::predict(plain, nd, type = "link")))
})
