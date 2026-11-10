# nolint start: object_usage_linter.
# ---- Fixtures ----
# A 20 by 20 grid of 30 m cells, all water, and a coverage table for two days
# whose p rises from west to east, with a strip of p = 0 in column 1.
avail_grid <- function() fixture_grid(20, 30)

avail_coverage <- function(grid = avail_grid()) {
  cells <- seq_len(terra::ncell(grid))
  col <- terra::colFromCell(grid, cells)
  p <- (col - 1) / 19
  rbind(data.frame(date = as.Date("2024-06-01"), cell = cells, p = p),
        data.frame(date = as.Date("2024-06-02"), cell = cells, p = rev(p)))
}

# Used points: animal a in two hours of 1 June, animal b in one hour of
# 2 June, three draws each, at cell centres.
avail_presence <- function(grid = avail_grid()) {
  t0 <- as.POSIXct("2024-06-01 00:00:00", tz = "UTC")
  cells <- c(105, 106, 125, 210, 211, 230, 300, 301, 320)
  xy <- terra::xyFromCell(grid, cells)
  data.frame(animal_id = rep(c("a", "a", "b"), each = 3),
             bin_start = t0 + rep(c(0, 3600, 86400 + 7200), each = 3),
             draw = rep(1:3, 3), cell = cells, x = xy[, 1], y = xy[, 2],
             stringsAsFactors = FALSE)
}

# ---- Draws and domains ----
test_that("each stratum draws from its own date's cells", {
  grid <- avail_grid()
  cover <- avail_coverage(grid)
  pres <- avail_presence(grid)
  av <- sample_available(pres, cover, grid, n = 500, seed = 1)
  # p rises west to east on 1 June and east to west on 2 June.
  col <- terra::colFromCell(grid, av$cell)
  june1 <- as.Date(av$bin_start, tz = "UTC") == as.Date("2024-06-01")
  expect_gt(mean(col[june1]), 12)
  expect_lt(mean(col[!june1]), 9)
})

test_that("cells with p = 0 are never drawn", {
  grid <- avail_grid()
  cover <- avail_coverage(grid)
  av <- sample_available(avail_presence(grid), cover, grid, n = 2000, seed = 2)
  p_of <- cover$p[match(paste(as.Date(av$bin_start, tz = "UTC"), av$cell),
                        paste(cover$date, cover$cell))]
  expect_true(all(p_of > 0))
  expect_false(any(terra::colFromCell(grid, av$cell[av$animal_id == "a"]) == 1))
})

test_that("weighted draws cover every cell of the date with p above 0", {
  grid <- avail_grid()
  cover <- avail_coverage(grid)
  pres <- avail_presence(grid)
  av <- sample_available(pres, cover, grid, n = 5000, seed = 3)
  for (d in as.Date(c("2024-06-01", "2024-06-02"))) {
    drawn <- unique(av$cell[as.Date(av$bin_start, tz = "UTC") == d])
    expect_setequal(drawn, cover$cell[cover$date == d & cover$p > 0])
  }
})

test_that("uniform draws take the cells at or above the threshold, evenly", {
  grid <- avail_grid()
  cover <- avail_coverage(grid)
  pres <- avail_presence(grid)[1, ]
  av <- sample_available(pres, cover, grid, n = 1e5, uniform_above = 0.5, seed = 4)
  day <- cover[cover$date == as.Date("2024-06-01"), ]
  eligible <- day$cell[day$p >= 0.5]
  expect_true(all(av$cell %in% eligible))
  observed <- table(factor(av$cell, levels = eligible))
  test <- stats::chisq.test(as.vector(observed))
  expect_gt(test$p.value, 0.001)
  expect_identical(attr(av, "uniform_above"), 0.5)
  expect_match(attr(av, "estimand"), "uniform over the cells with p at or above 0.5")
  # uniform_above = 0 takes the whole water body, p = 0 cells included.
  all_water <- sample_available(pres, cover, grid, n = 2e4, uniform_above = 0, seed = 4)
  expect_setequal(unique(all_water$cell), day$cell)
})

test_that("n per stratum is exact and the strata follow the presence draws", {
  grid <- avail_grid()
  pres <- avail_presence(grid)
  av <- sample_available(pres, avail_coverage(grid), grid, n = 7, seed = 4)
  per <- table(paste(av$animal_id, av$bin_start, av$draw))
  expect_equal(as.vector(per), rep(7L, 9))
  expect_setequal(names(per), paste(pres$animal_id, pres$bin_start, pres$draw))
  expect_named(av, c("animal_id", "bin_start", "draw", "cell", "x", "y"))
  expect_equal(cbind(av$x, av$y), unname(terra::xyFromCell(grid, av$cell)))
  expect_identical(attr(av, "n"), 7L)
  expect_null(attr(av, "uniform_above"))
  expect_identical(attr(av, "epsg"), 32617L)
  expect_match(attr(av, "estimand"), "weighted by the probability that the bin is detected")
  expect_identical(nrow(attr(av, "dropped")), 0L)
})

test_that("weighted draw frequencies follow p", {
  grid <- avail_grid()
  cover <- avail_coverage(grid)
  pres <- avail_presence(grid)[1, ]
  av <- sample_available(pres, cover, grid, n = 1e5, seed = 5)
  day <- cover[cover$date == as.Date("2024-06-01") & cover$p > 0, ]
  expected <- day$p / sum(day$p)
  observed <- table(factor(av$cell, levels = day$cell))
  expect_identical(sum(observed), 1e5L)
  test <- suppressWarnings(stats::chisq.test(as.vector(observed), p = expected))
  expect_gt(test$p.value, 0.001)
})

test_that("a seed gives the same draws and leaves the caller's stream alone", {
  grid <- avail_grid()
  pres <- avail_presence(grid)
  cover <- avail_coverage(grid)
  set.seed(99)
  before <- .Random.seed
  a <- sample_available(pres, cover, grid, n = 10, seed = 5)
  expect_identical(.Random.seed, before)
  expect_identical(sample_available(pres, cover, grid, n = 10, seed = 5), a)
  expect_false(identical(sample_available(pres, cover, grid, n = 10, seed = 6),
                         a))
  rm(".Random.seed", envir = globalenv())
  sample_available(pres, cover, grid, n = 1, seed = 5)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_s3_class(sample_available(pres, cover, grid, n = 1), "data.frame")
})

test_that("empty strata are dropped, counted and warned", {
  grid <- avail_grid()
  pres <- avail_presence(grid)
  cover <- avail_coverage(grid)
  # No coverage on 2 June: animal b's three strata have no date.
  one_day <- cover[cover$date == as.Date("2024-06-01"), ]
  expect_warning(av <- sample_available(pres, one_day, grid, n = 4, seed = 1),
                 "3 of 9 strata .* 3 with no coverage on their date, 0 with no cell")
  expect_identical(nrow(av), 24L)
  dropped <- attr(av, "dropped")
  expect_identical(dropped$animal_id, rep("b", 3))
  expect_identical(dropped$reason, rep("no coverage on the date", 3))
  # Every cell of 2 June at p = 0: animal b's strata have no qualifying cell.
  blank <- cover
  blank$p[blank$date == as.Date("2024-06-02")] <- 0
  expect_warning(av2 <- sample_available(pres, blank, grid, n = 4, seed = 1),
                 "0 with no coverage on their date, 3 with no cell of p\\s+above\\s+0")
  expect_identical(unique(attr(av2, "dropped")$reason), "no qualifying cell on the date")
  # Halve p: no cell reaches 0.9, so every stratum goes.
  expect_warning(sample_available(pres, transform(cover, p = p / 2), grid, n = 4,
                                  uniform_above = 0.9, seed = 1),
                 "9 of 9 strata.*no\\s+cell\\s+of\\s+p\\s+at\\s+or\\s+above\\s+0.9")
})

test_that("draws from a location surface work with system_de()", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80,
                         key = "receiver_serial")
  counts <- rbind(counts_for(field, c(3, 1, 0, 0), 12, "a", 0),
                  counts_for(field, c(0, 4, 0, 0), 12, "b", 1))
  presence <- sample_locations(location_surface(counts, field), 20, seed = 1)
  listening <- data.frame(receiver_serial = field$receivers, date = as.Date("2024-06-01"))
  cover <- system_de(field, listening, transmissions = 24)
  av <- sample_available(presence, cover, field$grid, n = 5, seed = 1)
  expect_identical(nrow(av), nrow(presence) * 5L)
  expect_true(all(av$cell %in% cover$cell[cover$p > 0]))
})

# ---- Checks ----
test_that("bad arguments stop", {
  grid <- avail_grid()
  pres <- avail_presence(grid)
  cover <- avail_coverage(grid)
  sa <- function(p = pres, cv = cover, g = grid, n = 2, uniform_above = NULL) {
    sample_available(p, cv, g, n = n, uniform_above = uniform_above)
  }
  expect_error(sa(n = 0), "`n` must be one whole number")
  expect_error(sa(n = 2.5), "`n` must be one whole number")
  expect_error(sa(uniform_above = 1.5), "`uniform_above` must be NULL or one number in 0 to 1")
  expect_error(sa(uniform_above = NA), "`uniform_above` must be NULL or one number in 0 to 1")
  expect_error(sa(uniform_above = c(0.1, 0.2)), "`uniform_above` must be NULL")
  expect_error(sa(g = "grid"), "`grid` must be a SpatRaster")
  expect_error(sa(p = pres[0, ]), "`presence` has no rows")
  expect_error(sa(p = pres[, -2]), "`presence` has no column `bin_start`")
  expect_error(sa(p = transform(pres, animal_id = factor(animal_id))),
               "Column `animal_id` of `presence` must be character, not factor")
  local_pres <- pres
  attr(local_pres$bin_start, "tzone") <- "America/Toronto"
  expect_error(sa(p = local_pres), "must be in UTC")
  expect_error(sa(p = transform(pres, draw = 0)), "whole numbers of at least 1")
  expect_error(sa(p = transform(pres, x = NA_real_)), "must be finite")
  expect_error(sa(p = rbind(pres, pres[1, ])), "1 row repeats")
  expect_error(sa(p = transform(pres, cell = cell + 1)), "9 of 9 rows of `presence` have a")
  epsg_pres <- pres
  attr(epsg_pres, "epsg") <- 4326L
  expect_error(sa(p = epsg_pres), "EPSG:4326 and the grid in EPSG:32617")
  expect_error(sa(cv = cover[0, ]), "`coverage` has no rows")
  expect_error(sa(cv = cover[, -3]), "`coverage` has no column `p`")
  expect_error(sa(cv = transform(cover, p = p * 2)), "must lie in 0 to 1")
  expect_error(sa(cv = transform(cover, cell = cell + 1000)), "cells of `grid`, 1 to 400")
  expect_error(sa(cv = rbind(cover, cover[1, ])), "one row per date and cell")
  expect_error(sa(cv = transform(cover, date = format(date))), "must be Date")
})

# nolint end
