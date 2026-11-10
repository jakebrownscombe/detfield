test_that("the path around the peninsula is the closed-form shortest path", {
  grid <- peninsula_lake()
  d <- water_distances(grid, peninsula_receivers())
  # From (285, 435) to (885, 435) the path bends round the peninsula tip at
  # the corners (570, 180) and (600, 180): 2 * sqrt(285^2 + 255^2) + 30 m.
  across <- pair(d, "west", terra::cellFromRowCol(grid, 6, 30))
  expect_equal(across$distance_m, 2 * sqrt(285^2 + 255^2) + 30, tolerance = 1e-9)
  expect_equal(across$straight_m, 600)
  expect_false(across$in_water)
  # South of the tip the straight line still crosses the peninsula.
  south <- pair(d, "west", terra::cellFromRowCol(grid, 20, 30))
  expect_false(south$in_water)
  expect_gt(south$distance_m, south$straight_m)
  # On the same side the straight line is used.
  near <- pair(d, "west", terra::cellFromRowCol(grid, 18, 2))
  expect_true(near$in_water)
  expect_equal(near$distance_m, sqrt((45 - 285)^2 + (75 - 435)^2))
  # The grid path round the peninsula, from a hand count on 30 m cells.
  g <- water_distances(grid, peninsula_receivers(), control = list(path = "grid"))
  expect_equal(pair(g, "west", terra::cellFromRowCol(grid, 6, 30))$distance_m,
               2 * (9 * sqrt(2) + 1) * 30, tolerance = 1e-9)
})

test_that("the table has one row per receiver and water cell, and its attributes", {
  grid <- peninsula_lake()
  d <- water_distances(grid, peninsula_receivers())
  n_water <- sum(!is.na(terra::values(grid)))
  expect_s3_class(d, c("water_distances", "data.frame"))
  expect_identical(names(d), c("cell", "station_id", "distance_m", "straight_m", "in_water"))
  expect_equal(nrow(d), 2 * n_water)
  expect_equal(attr(d, "cells"), which(!is.na(terra::values(grid, mat = FALSE))))
  expect_identical(attr(d, "key"), "station_id")
  expect_null(attr(d, "lookup"))
  expect_equal(attr(d, "unreachable"), 0)
  expect_true(all(d$distance_m >= d$straight_m - 1e-9 | d$in_water))
  expect_equal(d$distance_m[d$in_water], d$straight_m[d$in_water])
  expect_s4_class(attr(d, "grid"), "PackedSpatRaster")
})

test_that("max_distance drops far pairs and unreachable water is counted", {
  grid <- peninsula_lake()
  # A pond cut off by land: cells (18:20, 36:40) ringed by a wall at column 35
  # and row 17.
  grid[terra::cellFromRowColCombine(grid, 17:20, 35)] <- NA
  grid[terra::cellFromRowColCombine(grid, 17, 35:40)] <- NA
  d <- water_distances(grid, peninsula_receivers()[1, ])
  pond <- terra::cellFromRowColCombine(grid, 18:20, 36:40)
  expect_false(any(pond %in% d$cell))
  expect_equal(attr(d, "unreachable"), length(pond))
  near <- water_distances(grid, peninsula_receivers()[1, ], max_distance = 300)
  expect_lte(max(near$distance_m), 300)
  expect_equal(nrow(near), sum(d$distance_m <= 300))
  expect_error(water_distances(grid, peninsula_receivers(), max_distance = -1), "positive")
})

test_that("the line test is exact at corners", {
  grid <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 90, ymin = 0, ymax = 90,
                      crs = "EPSG:32617", vals = 1)
  grid[terra::cellFromRowCol(grid, 2, 2)] <- NA
  info <- grid_info(grid)
  # Through the centre of the land cell: blocked.
  expect_false(line_clear(info, 15, 45, 75, 45))
  # The land cell spans 30 to 60 m on both axes. The line x + y = 120 touches
  # only its corner at (60, 60), so it stays clear; 1 m lower it enters.
  expect_true(line_clear(info, 45, 75, 75, 45))
  expect_false(line_clear(info, 45, 74, 75, 44))
  expect_false(line_clear(info, 15, 15, 75, 75))
  expect_true(line_clear(info, 15, 15, 15, 15))
})

test_that("a line through the point where two land cells meet diagonally is blocked", {
  # 30 m cells; land at row 2, column 2 (x 30 to 60, y 60 to 90) and row 3,
  # column 3 (x 60 to 90, y 30 to 60). They meet at the corner (60, 60).
  grid <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 120, ymin = 0, ymax = 120,
                      crs = "EPSG:32617", vals = 1)
  grid[terra::cellFromRowColCombine(grid, 2, 2)] <- NA
  grid[terra::cellFromRowColCombine(grid, 3, 3)] <- NA
  info <- grid_info(grid)
  # The line from the centre at (75, 75) to the centre at (15, 15) passes the
  # corner between the two land cells: the gap is closed.
  expect_false(line_clear(info, 75, 75, 15, 15))
  expect_false(line_clear(info, 15, 15, 75, 75))
  expect_false(line_clear(info, 105, 105, 15, 15))
  # Shifted a few millimetres either way it passes through one land cell.
  for (shift in c(-0.004, 0.004)) {
    expect_false(line_clear(info, 75 + shift, 75, 15, 15))
    expect_false(line_clear(info, 75, 75 + shift, 15, 15))
  }
  # Along the other diagonal it crosses both land cells' centres.
  expect_false(line_clear(info, 15, 105, 105, 15))
  # Past the corner by more than the tolerance (0.03 m on 30 m cells) the
  # line is decided by the cells it enters: 1 m off, it enters land.
  expect_false(line_clear(info, 76, 75, 16, 15))
})

test_that("a line that grazes one land corner is clear however the receiver is rounded", {
  grid <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 90, ymin = 0, ymax = 90,
                      crs = "EPSG:32617", vals = 1)
  grid[terra::cellFromRowCol(grid, 2, 2)] <- NA
  info <- grid_info(grid)
  # The line x + y = 120 touches only the corner (60, 60) of the land cell.
  shifts <- expand.grid(dx = c(-0.006, -0.002, 0, 0.002, 0.006), dy = c(-0.006, 0, 0.006))
  for (i in seq_len(nrow(shifts))) {
    expect_true(line_clear(info, 45 + shifts$dx[i], 75 + shifts$dy[i], 75, 45))
    expect_true(line_clear(info, 75 + shifts$dx[i], 45 + shifts$dy[i], 45, 75))
  }
  # The tolerance is a thousandth of the cell side.
  expect_equal(touch_tolerance(info), 0.03)
  expect_true(line_clear(info, 45, 75 - 0.029, 75, 45 - 0.029))
  expect_false(line_clear(info, 45, 75 - 0.5, 75, 45 - 0.5))
})

test_that("a line along a grid line is blocked only where land lies on both sides", {
  # Land west of x = 30 in row 2, and east of it in rows 2 and 3.
  grid <- terra::rast(nrows = 4, ncols = 2, xmin = 0, xmax = 60, ymin = 0, ymax = 120,
                      crs = "EPSG:32617", vals = 1)
  info <- grid_info(grid)
  grid_one_side <- grid
  grid_one_side[terra::cellFromRowColCombine(grid, 2:3, 2)] <- NA
  expect_true(line_clear(grid_info(grid_one_side), 30, 105, 30, 15))
  grid_both <- grid_one_side
  grid_both[terra::cellFromRowColCombine(grid, 2, 1)] <- NA
  expect_false(line_clear(grid_info(grid_both), 30, 105, 30, 15))
  expect_true(line_clear(info, 30, 105, 30, 15))
  # Horizontal: land north of y = 60 in column 1 and south of it in column 2.
  grid_h <- terra::rast(nrows = 2, ncols = 4, xmin = 0, xmax = 120, ymin = 0, ymax = 60,
                        crs = "EPSG:32617", vals = 1)
  grid_h[terra::cellFromRowColCombine(grid_h, 1, 2)] <- NA
  expect_true(line_clear(grid_info(grid_h), 15, 30, 105, 30))
  grid_h[terra::cellFromRowColCombine(grid_h, 2, 3)] <- NA
  expect_false(line_clear(grid_info(grid_h), 15, 30, 105, 30))
})

test_that("a receiver on a grid line is tested on the side its line leaves by", {
  # Land at row 2, column 1 (x 0 to 30, y 30 to 60). The receiver sits on the
  # grid line x = 30, in column 2 by floor(), and its line to the centre at
  # (15, 75) runs through the land cell from its first metre.
  grid <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 90, ymin = 0, ymax = 90,
                      crs = "EPSG:32617", vals = 1)
  grid[terra::cellFromRowCol(grid, 2, 1)] <- NA
  info <- grid_info(grid)
  expect_false(line_clear(info, 30, 45, 15, 75))
  expect_false(line_clear(info, 30.004, 45, 15, 75))
  expect_false(line_clear(info, 29.996, 45, 15, 75))
  # Leaving east, the land cell behind the receiver is not on the line.
  expect_true(line_clear(info, 30, 45, 75, 75))
})

test_that("receivers a few millimetres off their cell centres give the same blocked pairs and field", { # nolint: line_length_linter.
  # Station files store latitude and longitude to seven decimals, which moves
  # a station up to about 1 cm off the centre it was placed at.
  set.seed(7)
  grid <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 1200, ymin = 0, ymax = 1200,
                      crs = "EPSG:32617", vals = 1)
  grid[sample(terra::ncell(grid), 240)] <- NA
  info <- grid_info(grid)
  at <- xy_of_cell(info, sample(which(info$water), 12))
  centred <- data.frame(station_id = sprintf("S%02d", 1:12), x = at[, 1], y = at[, 2],
                        deploy_z = 5)
  model <- de_logistic(d50 = 250, d05 = 500)
  base <- water_distances(grid, centred)
  base_field <- suppressMessages(de_field(base, model))
  key <- function(d) sort(paste(d$station_id, d$cell)[!d$in_water])
  for (dir in list(c(1, 0), c(-1, 0), c(0, 1), c(0, -1), c(1, 1), c(-1, 1), c(1, -1),
                   c(-1, -1))) {
    moved <- centred
    moved$x <- moved$x + 0.006 * dir[1]
    moved$y <- moved$y + 0.004 * dir[2]
    d <- water_distances(grid, moved)
    expect_identical(key(d), key(base))
    f <- suppressMessages(de_field(d, model))
    expect_equal(f$p, base_field$p, tolerance = 1e-4)
  }
  # Each station moved its own way, as rounding does.
  moved <- centred
  moved$x <- moved$x + stats::runif(12, -0.007, 0.007)
  moved$y <- moved$y + stats::runif(12, -0.007, 0.007)
  expect_identical(key(water_distances(grid, moved)), key(base))
})

test_that("Stoney receivers a few millimetres off their centres keep their blocked pairs", {
  grid <- terra::rast(system.file("extdata", "stoney_depth.tif", package = "detfield"))
  rx <- stoney_track$receivers[seq(1, 70, by = 7), c("station_id", "x", "y", "deploy_z")]
  base <- water_distances(grid, rx)
  key <- function(d) sort(paste(d$station_id, d$cell)[!d$in_water])
  set.seed(11)
  moved <- rx
  moved$x <- moved$x + stats::runif(nrow(rx), -0.0065, 0.0065)
  moved$y <- moved$y + stats::runif(nrow(rx), -0.0065, 0.0065)
  d <- water_distances(grid, moved)
  expect_identical(key(d), key(base))
  expect_lt(max(abs(d$distance_m - base$distance_m)), 0.01)
})

test_that("the line test agrees with dense sampling on a random lake", {
  set.seed(42)
  grid <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 1200, ymin = 0, ymax = 1200,
                      crs = "EPSG:32617", vals = 1)
  grid[sample(terra::ncell(grid), 160)] <- NA
  info <- grid_info(grid)
  water <- which(info$water)
  from <- xy_of_cell(info, sample(water, 1))
  to <- xy_of_cell(info, sample(water, 400))
  exact <- line_clear(info, from[1], from[2], to[, 1], to[, 2])
  dense <- vapply(seq_len(nrow(to)), function(i) {
    t <- seq(0, 1, length.out = 20000)
    cells <- cell_of_xy(info, from[1] + t * (to[i, 1] - from[1]),
                        from[2] + t * (to[i, 2] - from[2]))
    all(info$water[cells])
  }, logical(1))
  expect_identical(exact, dense)
})

test_that("receivers come as x and y, as latitude and longitude, or as sf", {
  grid <- peninsula_lake()
  rx <- peninsula_receivers()
  ll <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(rx, coords = c("x", "y"), crs = 32617),
                                            4326))
  rx_ll <- data.frame(station_id = rx$station_id, deploy_lat = ll[, 2], deploy_lon = ll[, 1])
  a <- water_distances(grid, rx)
  b <- water_distances(grid, rx_ll)
  expect_equal(b$distance_m, a$distance_m, tolerance = 1e-6)
  c <- water_distances(grid, sf::st_as_sf(rx, coords = c("x", "y"), crs = 32617))
  expect_equal(c$distance_m, a$distance_m)
  attr(rx, "epsg") <- 32618L
  expect_error(water_distances(grid, rx), "epsg 32618")
  attr(rx, "epsg") <- 32617L
  expect_silent(water_distances(grid, rx))
})

test_that("the grid can be a file path and must be projected in metres", {
  grid <- peninsula_lake()
  path <- withr::local_tempfile(fileext = ".tif")
  terra::writeRaster(grid, path)
  expect_equal(nrow(water_distances(path, peninsula_receivers())),
               nrow(water_distances(grid, peninsula_receivers())))
  lonlat <- terra::rast(nrows = 5, ncols = 5, crs = "EPSG:4326", vals = 1)
  err <- expect_error(water_distances(lonlat, peninsula_receivers()), "longitude and latitude",
                      class = "detfield_error_grid")
  expect_s3_class(err, "detfield_error")
  expect_identical(conditionCall(err), quote(water_distances(lonlat, peninsula_receivers())))
  two <- c(grid, grid * 2)
  expect_equal(nrow(water_distances(two, peninsula_receivers())),
               nrow(water_distances(grid, peninsula_receivers())))
})

test_that("each grid refusal carries the grid class and the user's call", {
  rx <- peninsula_receivers()
  feet_ny <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 500, ymin = 0, ymax = 500,
                         crs = "EPSG:2263", vals = 1)
  feet_ca <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 500, ymin = 0, ymax = 500,
                         crs = "EPSG:2227", vals = 1)
  km <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5,
                    crs = "+proj=utm +zone=17 +datum=WGS84 +units=km", vals = 1)
  no_crs <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 500, ymin = 0, ymax = 500,
                        crs = "", vals = 1)
  lonlat <- terra::rast(nrows = 5, ncols = 5, crs = "EPSG:4326", vals = 1)
  cases <- list(
    list(call = quote(water_distances(list(1), rx)), text = "must be a SpatRaster"),
    list(call = quote(water_distances(no_crs, rx)), text = "has no CRS"),
    list(call = quote(water_distances(lonlat, rx)), text = "longitude and latitude"),
    list(call = quote(water_distances(feet_ny, rx)), text = "0.3048006 m"),
    list(call = quote(water_distances(feet_ca, rx)), text = "0.3048006 m"),
    list(call = quote(water_distances(km, rx)), text = "is 1000 m")
  )
  for (case in cases) {
    err <- expect_error(eval(case$call), case$text, fixed = TRUE,
                        class = "detfield_error_grid")
    expect_s3_class(err, "detfield_error")
    expect_identical(conditionCall(err), case$call)
  }
  caught <- tryCatch(water_distances(feet_ny, rx), detfield_error_grid = function(e) "caught")
  expect_identical(caught, "caught")
})

test_that("a CRS whose scale at the raster centre is far from 1 is refused", {
  near <- function(lon, lat, crs) {
    c0 <- terra::project(cbind(lon, lat), from = "EPSG:4326", to = crs)
    terra::rast(nrows = 10, ncols = 10, xmin = c0[1] - 500, xmax = c0[1] + 500,
                ymin = c0[2] - 500, ymax = c0[2] + 500, crs = crs, vals = 1)
  }
  rx <- peninsula_receivers()
  # Web Mercator at 45 N: about 1 / cos(45) = 1.414 on each segment; its
  # spherical formulas on an ellipsoid make the two axes differ in the third
  # decimal.
  merc <- near(-81, 45, "EPSG:3857")
  s <- crs_scale(merc)
  expect_equal(s$scale, rep(1 / cos(45 * pi / 180), 4), tolerance = 0.01)
  err <- expect_error(water_distances(merc, rx), class = "detfield_error_grid")
  expect_s3_class(err, "detfield_error")
  expect_identical(conditionCall(err), quote(water_distances(merc, rx)))
  expect_match(conditionMessage(err), "1.41[0-9] east-west,[[:space:]]+1.41[0-9] north-south")
  expect_match(conditionMessage(err), "17N (EPSG:32617)", fixed = TRUE)
  # Sinusoidal off its central meridian: both axes within 1 percent, the
  # diagonals 0.934 and 1.071.
  sinus <- near(11, 45.6, "ESRI:54008")
  s <- crs_scale(sinus)
  expect_true(all(abs(s$scale[1:2] - 1) <= max_scale_error))
  expect_equal(s$scale[3:4], c(0.934, 1.071), tolerance = 0.001)
  err <- expect_error(water_distances(sinus, rx), class = "detfield_error_grid")
  expect_match(conditionMessage(err), "0.934 south-west to north-east")
  expect_match(conditionMessage(err), "32N (EPSG:32632)", fixed = TRUE)
  # The same CRS on its central meridian, and UTM in its zone, pass.
  expect_s4_class(as_water_grid(near(0, 45.6, "ESRI:54008")), "SpatRaster")
  expect_s4_class(as_water_grid(near(-81, 45, "EPSG:32617")), "SpatRaster")
  # South of the equator the zone takes the 327zz code.
  err <- expect_error(as_water_grid(near(150, -33, "EPSG:3857")), class = "detfield_error_grid")
  expect_match(conditionMessage(err), "56S (EPSG:32756)", fixed = TRUE)
  # The unit refusal comes before the scale is measured.
  feet <- near(-74, 40.7, "EPSG:2263")
  expect_error(as_water_grid(feet), "0.3048006 m", fixed = TRUE, class = "detfield_error_grid")
})

test_that("a CRS whose scale cannot be measured is refused with no other condition", {
  rx <- peninsula_receivers()
  lake <- peninsula_lake()
  eng <- lake
  terra::crs(eng) <- paste0('ENGCRS["Lake grid",EDATUM["Arbitrary"],CS[Cartesian,2],',
                            'AXIS["easting",east,LENGTHUNIT["metre",1]],',
                            'AXIS["northing",north,LENGTHUNIT["metre",1]]]')
  expect_true(all(is.na(crs_scale(eng)$scale)))
  warnings <- 0L
  err <- withCallingHandlers(
    tryCatch(water_distances(eng, rx), error = identity),
    warning = function(w) {
      warnings <<- warnings + 1L
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(warnings, 0L)
  expect_s3_class(err, "detfield_error_grid")
  expect_s3_class(err, "detfield_error")
  expect_identical(conditionCall(err), quote(water_distances(eng, rx)))
  expect_match(conditionMessage(err), "cannot be measured")
  expect_match(conditionMessage(err), "projected CRS with a geographic datum")
  # The same grid in its own projected CRS still passes.
  expect_s3_class(water_distances(lake, rx), "water_distances")
})

test_that("bad receiver tables stop with the column at fault", {
  grid <- peninsula_lake()
  rx <- peninsula_receivers()
  expect_error(water_distances(grid, rx[, c("x", "y")]), "station_id or a receiver_serial")
  bad <- rx
  bad$station_id <- factor(bad$station_id)
  expect_error(water_distances(grid, bad), "station_id` of `receivers` must be character")
  bad <- rx
  bad$x[1] <- NA
  expect_error(water_distances(grid, bad), "Column `x`")
  expect_error(water_distances(grid, rx[, c("station_id", "deploy_z")]), "deploy_lat")
  expect_error(water_distances(grid, as.list(rx)), "data frame or an sf")
  expect_error(water_distances(grid, rx, control = list(cores = 1, speed = 2)), "speed")
  expect_error(water_distances(grid, rx, control = 1), "must be a list")
  expect_error(water_distances(grid, rx, control = list(path = "lattice")),
               "`control\\$path` must be")
})

test_that("deployments of one station collapse to their mean position", {
  grid <- peninsula_lake()
  rx <- data.frame(station_id = c("west", "west", "east"), receiver_serial = c("A", "B", "A"),
                   x = c(280, 290, 885), y = c(435, 435, 435),
                   deploy_datetime = as.POSIXct(c("2024-01-01", "2024-06-01", "2024-06-01"),
                                                tz = "UTC"),
                   recover_datetime = as.POSIXct(c("2024-06-01", NA, NA), tz = "UTC"))
  d <- water_distances(grid, rx)
  r <- attr(d, "receivers")
  expect_equal(r$x[r$station_id == "west"], 285)
  expect_equal(r$n_rows, c(2L, 1L))
  expect_equal(nrow(attr(d, "lookup")), 3)
  rx$x[2] <- 400
  expect_warning(water_distances(grid, rx), "more than one cell")
  rx$x[2] <- 290
  attr(rx$deploy_datetime, "tzone") <- "America/Toronto"
  expect_error(water_distances(grid, rx), "deploy_datetime` of `receivers` must be POSIXct in UTC")
})

test_that("a receiver on a land cell moves to water, or stops when far from it", {
  grid <- peninsula_lake()
  rx <- data.frame(receiver_serial = c("R1", "R2"), x = c(590, 885), y = c(500, 435))
  expect_warning(d <- water_distances(grid, rx), "nearest water cell, by .*: R1 \\(")
  r <- attr(d, "receivers")
  expect_identical(attr(d, "key"), "receiver_serial")
  expect_gt(r$snap_m[1], 0)
  expect_lte(r$snap_m[1], 30 * sqrt(2))
  expect_true(!is.na(terra::values(grid)[r$cell[1]]))
  far <- terra::rast(nrows = 20, ncols = 40, xmin = 0, xmax = 1200, ymin = 0, ymax = 600,
                     crs = "EPSG:32617", vals = 1)
  far[terra::cellFromRowColCombine(far, 1:10, 1:20)] <- NA
  expect_error(water_distances(far, data.frame(station_id = "a", x = 10, y = 590)),
               "from the nearest water cell")
  expect_error(water_distances(grid, data.frame(station_id = "a", x = -10, y = 5)),
               "outside `grid`")
})

test_that("several cores give the same table as one", {
  skip_on_os("windows")
  grid <- peninsula_lake()
  one <- water_distances(grid, peninsula_receivers())
  two <- water_distances(grid, peninsula_receivers(), control = list(cores = 2))
  expect_identical(as.data.frame(unclass(two)), as.data.frame(unclass(one)))
})

test_that("printing summarises the table", {
  d <- water_distances(peninsula_lake(), peninsula_receivers())
  expect_output(print(d), "1,572 pairs: 2 receivers \\(by station_id\\) x 786 water cells",
                width = 200)
  expect_output(print(d), "Straight line used for")
  plain <- d
  attr(plain, "receivers") <- NULL
  expect_output(print(plain), "distance_m")
})

test_that("an error in a worker process reaches the caller", {
  skip_on_os("windows")
  fail_on_two <- function(i) if (i == 2) stop("worker ", i, " failed") else i
  expect_error(parallel_apply(1:2, fail_on_two, cores = 2), "worker 2 failed")
  expect_identical(parallel_apply(1:2, identity, cores = 1), list(1L, 2L))
  classed <- function(i) if (i == 2) detfield_abort("bin 2", "degenerate") else i
  expect_error(parallel_apply(1:2, classed, cores = 2), "bin 2",
               class = "detfield_error_degenerate")
  expect_error(parallel_apply(1:2, function(i) if (i == 1) 1 else NULL, cores = 2),
               "ran out of memory", class = "detfield_error_worker")
})

test_that("the grid path stays within the documented 8.24% of the shortest water path", {
  # The 8-neighbour path is not the Euclidean
  # shortest path. Past the peninsula tip corners (570, 180) and (600, 180)
  # the shortest path is 2 * sqrt(285^2 + 255^2) + 30 = 794.85 m, which the
  # exact path gives.
  euclid <- 2 * sqrt(285^2 + 255^2) + 30
  cell <- terra::cellFromRowCol(peninsula_lake(), 6, 30)
  g <- water_distances(peninsula_lake(), peninsula_receivers(), control = list(path = "grid"))
  expect_gt(pair(g, "west", cell)$distance_m, euclid)
  expect_lte(pair(g, "west", cell)$distance_m / euclid, 1.0824)
  d <- water_distances(peninsula_lake(), peninsula_receivers())
  expect_equal(pair(d, "west", cell)$distance_m, euclid, tolerance = 1e-9)
})

test_that("distance_m is never shorter than the straight line", {
  # The grid path runs from the cell centre and the straight line from the
  # receiver, so near its cell's corner a receiver has 6 pairs whose grid
  # path is up to 4.25 m short; they are raised to the straight line. The
  # exact path starts at the receiver and needs no such step.
  rx <- data.frame(station_id = "a", x = 271, y = 449)
  g <- water_distances(peninsula_lake(), rx, control = list(path = "grid"))
  expect_true(all(g$distance_m >= g$straight_m))
  expect_equal(sum(!g$in_water & g$distance_m == g$straight_m), 6)
  d <- water_distances(peninsula_lake(), rx)
  expect_true(all(d$distance_m >= d$straight_m - 1e-9))
  expect_equal(sum(!d$in_water & d$distance_m == d$straight_m), 0)
})

test_that("land coded as a value rather than NA gives a warning", {
  # Land coded 0, or a grid with no NA cells, would otherwise be read as water.
  grid <- peninsula_lake()
  zero <- terra::classify(grid, cbind(NA, 0))
  expect_warning(water_distances(zero, peninsula_receivers()), "no NA cells, so all 800 cells")
  mixed <- grid
  mixed[terra::cellFromRowColCombine(mixed, 20, 1:40)] <- -5
  expect_warning(water_distances(mixed, peninsula_receivers()),
                 "40 of 786 water cells hold values of 0 or less and 746 hold values above 0")
  below <- grid * -1
  expect_silent(water_distances(below, peninsula_receivers()))
})

test_that("a grid with no CRS stops with a plain message", {
  # Without the check, the call stops with 'missing value where TRUE/FALSE needed'.
  grid <- peninsula_lake()
  terra::crs(grid) <- ""
  expect_error(water_distances(grid, peninsula_receivers()), "`grid` has no CRS")
})

test_that("any point of a land cell that borders water snaps", {
  # The limit is measured to the water edge, so a receiver 29.9 m from the
  # edge and 44.9 m from the nearest centre moves.
  grid <- terra::rast(nrows = 20, ncols = 40, xmin = 0, xmax = 1200, ymin = 0, ymax = 600,
                      crs = "EPSG:32617", vals = 1)
  grid[terra::cellFromRowColCombine(grid, 1:10, 1:20)] <- NA
  rx <- data.frame(station_id = "shore", x = 135, y = 329.9)
  expect_warning(d <- water_distances(grid, rx), "nearest water cell, by .*: shore \\(45 m\\)")
  expect_equal(attr(d, "receivers")$snap_m, 44.9, tolerance = 1e-9)
  expect_error(water_distances(grid, data.frame(station_id = "inland", x = 135, y = 372.5)),
               "inland is 72 m from the nearest water cell")
})

test_that("a serial at two positions stops instead of being averaged", {
  # Averaged, serial 123 at x = 285 and x = 885 would become one receiver at 585.
  rx <- data.frame(receiver_serial = c("123", "123"), x = c(285, 885), y = c(435, 435))
  expect_error(water_distances(peninsula_lake(), rx),
               "1 receiver_serial value sits at positions more than one cell \\(30 m\\)")
  close <- data.frame(receiver_serial = c("123", "123"), x = c(280, 290), y = c(435, 435))
  expect_equal(attr(water_distances(peninsula_lake(), close), "receivers")$x, 285)
})

# ---- Empty inputs ----
test_that("no receivers or no water stop with the cause named", {
  expect_error(water_distances(peninsula_lake(), peninsula_receivers()[0, ]),
               "`receivers` has no rows")
  dry <- peninsula_lake()
  terra::values(dry) <- NA
  expect_error(suppressWarnings(water_distances(dry, peninsula_receivers())),
               "`grid` has no water cells")
})

test_that("a snap can reach one and a half cell diagonals and no farther", {
  # A 2 by 2 block of water among land on 100 m cells; the receiver sits on
  # the diagonal from the block's corner, just inside one diagonal of its edge.
  grid <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000, ymin = 0, ymax = 1000,
                      crs = "EPSG:32617")
  m <- matrix(NA_real_, 10, 10)
  m[5:6, 1:2] <- 5
  terra::values(grid) <- as.vector(t(m))
  step <- 141.4 / sqrt(2)
  rx <- data.frame(station_id = c("in", "corner"), x = c(50, 200 + step), y = c(450, 600 + step))
  expect_warning(d <- water_distances(grid, rx), "corner \\(212 m\\)")
  snap <- attr(d, "receivers")$snap_m[2]
  expect_equal(snap, sqrt((150 - rx$x[2])^2 + (550 - rx$y[2])^2), tolerance = 1e-9)
  expect_gt(snap, 200)
  expect_lte(snap, 1.5 * sqrt(2) * 100)
  rx$x[2] <- 301
  rx$y[2] <- 701
  expect_error(water_distances(grid, rx), "from the nearest water cell")
})
