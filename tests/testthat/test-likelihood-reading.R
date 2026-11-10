# The two readings of the curve: "binomial" uses detections and misses at
# every listening receiver; "presence_only" weighs each cell by p for each
# detection and ignores misses.

# An open square of 41 by 41 cells of 50 m, with a logistic curve through
# 0.5 at 250 m and 0.05 at 500 m.
reading_field <- function(rx) {
  grid <- fixture_grid(side = 41, cell = 50)
  curve <- de_logistic(d50 = 250, d05 = 500)
  fixture_field(grid, rx, function(d) curve$predict(as.vector(d)) + 0 * d, sigma = 0)
}

reading_surface <- function(field, y, n, likelihood) {
  location_surface(counts_for(field, y, n), field, overdispersion = "none",
                   likelihood = likelihood, control = list(kept_mass = 1))
}

# The whole surface as one value per water cell, 0 where no cell was kept.
as_full <- function(s, field) {
  v <- numeric(terra::ncell(field$grid))
  v[s$cell] <- s$prob
  v[field$cells]
}

mode_distance <- function(s, field, rx) {
  xy <- terra::xyFromCell(field$grid, s$cell[which.max(s$prob)])
  sqrt((xy[1] - rx$x)^2 + (xy[2] - rx$y)^2)
}

test_that("one receiver hearing y of n: presence-only mode at the receiver, binomial on the ring", {
  rx <- data.frame(receiver_serial = c("R1", "R2"), x = c(1025, 1925), y = c(1025, 1925))
  field <- reading_field(rx)
  det <- reading_surface(field, c(15, 0), 30, "presence_only")
  bin <- reading_surface(field, c(15, 0), 30, "binomial")
  # Under "presence_only" the mode is the receiver's own cell.
  expect_equal(mode_distance(det, field, rx)[1], 0)
  # Under "binomial" the mode lies where p = 15 / 30 = 0.5, at d50 = 250 m,
  # within one cell of 50 m.
  expect_lt(abs(mode_distance(bin, field, rx)[1] - 250), 50)
  expect_identical(attr(det, "likelihood"), "presence_only")
  expect_identical(attr(bin, "likelihood"), "binomial")
})

test_that("two receivers: the mode sits between them under both readings", {
  rx <- data.frame(receiver_serial = c("R1", "R2"), x = c(825, 1225), y = c(1025, 1025))
  field <- reading_field(rx)
  # Between: nearer the midpoint (1025, 1025) than either receiver. Under the
  # binomial the two rings of p = 12 / 30 cross off the line joining them.
  for (lik in c("presence_only", "binomial")) {
    s <- reading_surface(field, c(12, 12), 30, lik)
    xy <- terra::xyFromCell(field$grid, s$cell[which.max(s$prob)])
    to_mid <- sqrt((xy[1] - 1025)^2 + (xy[2] - 1025)^2)
    expect_lt(to_mid, min(mode_distance(s, field, rx)), label = lik)
    expect_gt(xy[1], 825)
    expect_lt(xy[1], 1225)
  }
})

test_that("the presence-only reading ignores silent receivers", {
  rx <- data.frame(receiver_serial = c("R1", "R2"), x = c(1025, 1425), y = c(1025, 1025))
  rx3 <- rbind(rx, data.frame(receiver_serial = "R3", x = 1025, y = 1325))
  one <- reading_surface(reading_field(rx), c(15, 0), 30, "presence_only")
  two <- reading_surface(reading_field(rx3), c(15, 0, 0), 30, "presence_only")
  expect_equal(as_full(two, reading_field(rx3)), as_full(one, reading_field(rx)))
  # The binomial reading does use the silence.
  b1 <- reading_surface(reading_field(rx), c(15, 0), 30, "binomial")
  b2 <- reading_surface(reading_field(rx3), c(15, 0, 0), 30, "binomial")
  expect_gt(max(abs(as_full(b2, reading_field(rx3)) - as_full(b1, reading_field(rx)))), 1e-4)
})

test_that("overdispersion under the presence-only reading: none silently, station_hour warns", {
  rx <- data.frame(receiver_serial = c("R1", "R2"), x = c(1025, 1925), y = c(1025, 1925))
  field <- reading_field(rx)
  counts <- counts_for(field, c(15, 0), 30)
  expect_silent(location_surface(counts, field, overdispersion = "none",
                                 likelihood = "presence_only"))
  expect_silent(location_surface(counts, field, likelihood = "presence_only"))
  expect_warning(s <- location_surface(counts, field, overdispersion = "station_hour",
                                       likelihood = "presence_only"),
                 "ignored under likelihood")
  out <- paste(utils::capture.output(print(s)), collapse = "\n")
  expect_match(out, "Likelihood: presence-only")
  expect_match(out, "Spread: not used under the presence-only reading")
  expect_true("misfit_p" %in% names(attr(s, "bins")))
  expect_error(location_surface(counts, field, likelihood = "poisson"), "should be one of")
})
