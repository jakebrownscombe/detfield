# nolint start: object_usage_linter.
# ---- Fixtures ----
# The grid, coverage and presence of test-sample_available.R, rebuilt here,
# and two covariate layers on the grid: depth rising with the cell number and
# a column index.
ua_grid <- function() fixture_grid(20, 30)

ua_coverage <- function(grid = ua_grid()) {
  cells <- seq_len(terra::ncell(grid))
  p <- (terra::colFromCell(grid, cells) - 1) / 19
  rbind(data.frame(date = as.Date("2024-06-01"), cell = cells, p = p),
        data.frame(date = as.Date("2024-06-02"), cell = cells, p = rev(p)))
}

ua_presence <- function(grid = ua_grid()) {
  t0 <- as.POSIXct("2024-06-01 00:00:00", tz = "UTC")
  cells <- c(105, 106, 125, 210, 211, 230, 300, 301, 320)
  xy <- terra::xyFromCell(grid, cells)
  data.frame(animal_id = rep(c("a", "a", "b"), each = 3),
             bin_start = t0 + rep(c(0, 3600, 86400 + 7200), each = 3),
             draw = rep(1:3, 3), cell = cells, x = xy[, 1], y = xy[, 2],
             stringsAsFactors = FALSE)
}

ua_layers_fixture <- function(grid = ua_grid()) {
  depth <- grid
  terra::values(depth) <- seq_len(terra::ncell(grid)) / 10
  column <- grid
  terra::values(column) <- terra::colFromCell(grid, seq_len(terra::ncell(grid)))
  layers <- c(depth, column)
  names(layers) <- c("depth_m", "column")
  layers
}

# ---- Stacking ----
test_that("each stratum holds one used point and n available points", {
  grid <- ua_grid()
  pres <- ua_presence(grid)
  av <- sample_available(pres, ua_coverage(grid), grid, n = 6, seed = 1)
  st <- used_available(pres, av, ua_layers_fixture(grid))
  expect_named(st, c("stratum", "case", "animal_id", "bin_start", "draw", "cell", "x", "y",
                     "depth_m", "column"))
  expect_identical(nrow(st), 9L * 7L)
  expect_identical(as.vector(tapply(st$case, st$stratum, sum)), rep(1L, 9))
  expect_identical(as.vector(table(st$stratum[st$case == 0L])), rep(6L, 9))
  expect_identical(sort(unique(st$stratum)), 1:9)
  # The used point leads each stratum, and a stratum is one animal, bin and draw.
  expect_true(all(st$case[!duplicated(st$stratum)] == 1L))
  per <- unique(st[c("stratum", "animal_id", "bin_start", "draw")])
  expect_identical(nrow(per), 9L)
  # Used rows are the presence points.
  used <- st[st$case == 1L, ]
  expect_equal(used$cell, pres$cell[order(pres$animal_id, pres$bin_start, pres$draw)])
  expect_identical(sum(attr(st, "dropped")$rows), 0L)
  expect_identical(attr(st, "n"), 6L)
  expect_identical(attr(st, "epsg"), 32617L)
  expect_identical(attr(st, "estimand"), attr(av, "estimand"))
})

test_that("covariates equal the raster at the cell", {
  grid <- ua_grid()
  pres <- ua_presence(grid)
  layers <- ua_layers_fixture(grid)
  av <- sample_available(pres, ua_coverage(grid), grid, n = 50, seed = 2)
  st <- used_available(pres, av, layers)
  expect_equal(st$depth_m, terra::values(layers[["depth_m"]], mat = FALSE)[st$cell])
  expect_equal(st$column, terra::values(layers[["column"]], mat = FALSE)[st$cell])
  expect_equal(st$cell, terra::cellFromXY(layers, cbind(st$x, st$y)))
  # A point inside a cell, not at its centre, is read at that cell.
  off <- pres[, c("animal_id", "bin_start", "draw", "x", "y")]
  off$x <- off$x + 9
  off$y <- off$y - 11
  st_off <- used_available(off, av, layers)
  expect_equal(st_off$depth_m[st_off$case == 1L], st$depth_m[st$case == 1L])
  # A categorical layer keeps its levels.
  habitat <- terra::classify(layers[["column"]], cbind(c(0, 10), c(10, 21), c(1, 2)))
  levels(habitat) <- data.frame(id = 1:2, habitat = c("west", "east"))
  names(habitat) <- "habitat"
  st_cat <- used_available(pres, av, habitat)
  expect_s3_class(st_cat$habitat, "factor")
  expect_setequal(levels(st_cat$habitat), c("west", "east"))
})

test_that("rows with a missing layer value are dropped and counted", {
  grid <- ua_grid()
  pres <- ua_presence(grid)
  av <- sample_available(pres, ua_coverage(grid), grid, n = 30, seed = 3)
  layers <- ua_layers_fixture(grid)
  # Blank the depth of one used cell and of column 15, where no used point lies.
  blank_cells <- c(pres$cell[1], terra::cellFromRowColCombine(grid, 1:20, 15))
  expect_false(any(terra::colFromCell(grid, pres$cell) == 15))
  layers[["depth_m"]][blank_cells] <- NA
  in_lost <- av$animal_id == "a" & av$draw == 1 & av$bin_start == pres$bin_start[1]
  n_avail_na <- sum(av$cell %in% blank_cells & !in_lost)
  n_lost <- sum(av$animal_id == "a" & av$draw == 1 & av$bin_start == pres$bin_start[1])
  expect_warning(st <- used_available(pres, av, layers), "1 used and .* rows dropped")
  dropped <- attr(st, "dropped")
  expect_identical(dropped$rows[dropped$case == 1L & dropped$reason == "missing layer value"],
                   1L)
  expect_identical(dropped$rows[dropped$case == 0L & dropped$reason == "missing layer value"],
                   n_avail_na)
  expect_identical(dropped$rows[dropped$reason == "used point has a missing layer value"],
                   n_lost)
  expect_false(anyNA(st$depth_m))
  expect_identical(nrow(st), 9L * 31L - 1L - n_avail_na - n_lost)
  expect_identical(length(unique(st$stratum)), 8L)
  expect_identical(as.vector(tapply(st$case, st$stratum, sum)), rep(1L, 8))
})

test_that("strata without available points are dropped and counted", {
  grid <- ua_grid()
  pres <- ua_presence(grid)
  cover <- ua_coverage(grid)
  expect_warning(av <- sample_available(pres, cover[cover$date == as.Date("2024-06-01"), ],
                                        grid, n = 4, seed = 1), "3 of 9 strata")
  expect_warning(st <- used_available(pres, av, ua_layers_fixture(grid)),
                 "3 used and 0 available rows dropped: 3 used points with no")
  expect_identical(attr(st, "dropped")$rows[1], 3L)
  expect_identical(nrow(st), 6L * 5L)
  # A stratum whose every available point lacks a value goes whole: move its
  # available points to cell 400, far from every other point, and blank it.
  layers <- ua_layers_fixture(grid)
  first <- av$animal_id == "a" & av$draw == 2 & av$bin_start == min(av$bin_start)
  expect_false(400 %in% c(av$cell, pres$cell))
  av$cell[first] <- 400
  far <- terra::xyFromCell(grid, 400)
  av$x[first] <- far[1, 1]
  av$y[first] <- far[1, 2]
  layers[["column"]][400] <- NA
  st2 <- suppressWarnings(used_available(pres, av, layers))
  d2 <- attr(st2, "dropped")
  expect_identical(d2$rows[d2$reason == "every available point has a missing layer value"], 1L)
  expect_identical(length(unique(st2$stratum)), 5L)
})

# ---- Checks ----
test_that("bad arguments stop", {
  grid <- ua_grid()
  pres <- ua_presence(grid)
  av <- sample_available(pres, ua_coverage(grid), grid, n = 2, seed = 1)
  layers <- ua_layers_fixture(grid)
  expect_error(used_available(pres, av, "depth"), "`layers` must be a SpatRaster")
  bad_names <- layers
  names(bad_names) <- c("depth_m", "depth_m")
  expect_error(used_available(pres, av, bad_names), "unique names")
  clash <- layers
  names(clash) <- c("depth_m", "case")
  expect_error(used_available(pres, av, clash), "a layer named `case`")
  expect_error(used_available(pres, av[, -4], layers), "`available` has no column `cell`")
  expect_error(used_available(pres, av[0, ], layers), "`available` has no rows")
  expect_error(used_available(pres, transform(av, cell = cell + 1), layers),
               "rows of `available` have a")
  other <- transform(av, animal_id = "z")
  expect_error(used_available(pres, other, layers), "belong to no row of `presence`")
  expect_error(used_available(pres[, -1], av, layers), "`presence` has no column `animal_id`")
  epsg_av <- av
  attr(epsg_av, "epsg") <- 3857L
  expect_error(used_available(pres, epsg_av, layers), "`available` is in EPSG:3857")
})

# nolint end
