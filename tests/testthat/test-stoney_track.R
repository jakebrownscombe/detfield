stoney_grid <- function() {
  terra::rast(system.file("extdata", "stoney_depth.tif", package = "detfield", mustWork = TRUE))
}

stoney_field <- function(sigma = 0) {
  truth <- stoney_track$truth
  model <- de_model(de_logistic(truth$d50, truth$d05, depth = truth$depth),
                    sigma = list(station_hour = sigma))
  suppressMessages(de_field(water_distances(stoney_grid(), stoney_track$receivers), model))
}

test_that("the bundled raster is the Stoney depth raster", {
  path <- system.file("extdata", "stoney_depth.tif", package = "detfield", mustWork = TRUE)
  expect_identical(unname(tools::md5sum(path)), "21ba8e62d4f1fe7050b155a8fa6c59b3")
  v <- terra::values(stoney_grid(), mat = FALSE)
  # ?stoney_track and inst/extdata/README.md state these counts.
  expect_identical(sum(!is.na(v)), 4818L)
  expect_true(all(v <= 0, na.rm = TRUE))
  expect_identical(sum(v == 0, na.rm = TRUE), 56L)
})

test_that("the Stoney study runs from grid to regions", {
  field <- stoney_field()
  counts <- bin_detections(stoney_track$detections, stoney_track$receivers, stoney_track$tags)
  # 212 detected hours of five animals by 70 receivers, as
  # data-raw/stoney_track.R builds them.
  expect_equal(nrow(stoney_track$receivers), 70L)
  expect_equal(nrow(counts), 212L * 70L)
  expect_true(all(counts$n == 24L))
  surface <- location_surface(counts, field)
  region <- credible_region(surface)
  expect_equal(nrow(attr(region, "summary")), 212L)
  draws <- sample_locations(surface, n = 10, seed = 1)
  expect_equal(nrow(draws), 2120L)
})

test_that("regions hold a still animal on the Stoney field", {
  # The hourly likelihood describes an animal that stays in one cell for the
  # hour. Simulate such hours from the field itself: 40 cells heard by at
  # least one receiver, 24 transmissions, each receiver detecting each with
  # the field's probability. The 95% regions should hold the cell in about
  # 95% of hours.
  field <- stoney_field()
  heard <- which(apply(field$p, 1, max) >= 0.2)
  set.seed(11)
  rows <- sample(heard, 40)
  # One bin of the shipped counts gives the rows of every listening receiver
  shipped <- bin_detections(stoney_track$detections, stoney_track$receivers, stoney_track$tags)
  template <- shipped[shipped$animal_id == shipped$animal_id[1] &
                        shipped$bin_start == shipped$bin_start[1], ]
  template$animal_id <- "still"
  start <- template$bin_start[1]
  p_rx <- field$p[, match(template$station_id, field$receivers)]
  counts <- do.call(rbind, lapply(seq_along(rows), function(i) {
    hour <- template
    hour$bin_start <- start + (i - 1) * 3600
    hour$bin_end <- hour$bin_start + 3600
    hour$y <- stats::rbinom(nrow(hour), 24, p_rx[rows[i], ])
    hour$n <- 24L
    hour
  }))
  counts <- counts[ave(counts$y, counts$bin_start, FUN = sum) > 0, ]
  region <- credible_region(location_surface(counts, field))
  hours <- unique(region$bin_start)
  truth <- field$cells[rows][match(hours, start + (seq_along(rows) - 1) * 3600)]
  inside <- vapply(seq_along(hours), function(h) {
    truth[h] %in% region$cell[region$bin_start == hours[h]]
  }, NA)
  # The animal is heard in all 40 hours; the bound leaves room for two fewer.
  expect_gte(length(hours), 38L)
  expect_gte(mean(inside), 0.85)
})
