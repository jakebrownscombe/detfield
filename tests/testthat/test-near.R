# ---- Distance from the detecting receivers ----
# With de_logistic(300, 1000) and the plain binomial, a low count at one
# receiver is best explained by a ring of cells at the distance where the
# curve gives that rate, so an hourly mode can sit hundreds of metres from the
# receivers that heard the animal. This synthetic lake builds that case.

near_lake <- function() {
  grid <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 3000, ymin = 0, ymax = 3000,
                      crs = "EPSG:32617", vals = 1)
  grid[1] <- NA # one land cell in a far corner, so the grid has a shore
  rx <- data.frame(station_id = c("mid", "sw", "ne"), x = c(1525, 325, 2725),
                   y = c(1525, 325, 2725))
  field <- suppressMessages(de_field(water_distances(grid, rx), de_logistic(300, 1000)))
  list(grid = grid, rx = rx, field = field)
}

test_that("a ring away from the receiver is counted, though misfit_p passes it", {
  lake <- near_lake()
  ring <- counts_for(lake$field, c(1, 0, 0), 60, key = "station_id")
  close <- counts_for(lake$field, c(50, 0, 0), 60, start = 1, key = "station_id")
  s <- location_surface(rbind(ring, close), lake$field, "none", control = list(kept_mass = 1))
  b <- attr(s, "bins")
  expect_true(all(c("mode_cell", "mode_dist_m", "near_receiver", "near_m", "near_mass") %in%
                    names(b)))
  expect_equal(b$near_receiver, c("mid", "mid"))
  # The field's d50: cell centres sit on a 50 m lattice about the receiver,
  # so the longest clear distance at p >= 0.5 is 300 m exactly.
  expect_equal(b$near_m, c(300, 300))
  # The ring: almost no mass near the receiver, a mode over 1 km away, and a
  # misfit tail that does not flag it.
  expect_lt(b$near_mass[1], 0.001)
  expect_gt(b$mode_dist_m[1], 1000)
  expect_gt(b$misfit_p[1], 0.1)
  # The strong count: most mass near the receiver.
  expect_gt(b$near_mass[2], 0.5)
  expect_lt(b$mode_dist_m[2], 300)
  # Independent recount from the surface rows and the grid.
  xy <- terra::xyFromCell(lake$grid, s$cell)
  d_mid <- sqrt((xy[, 1] - 1525)^2 + (xy[, 2] - 1525)^2)
  for (i in 1:2) {
    one <- as.numeric(s$bin_start) == as.numeric(b$bin_start[i])
    expect_equal(b$near_mass[i], sum(s$prob[one & d_mid <= 300]), tolerance = 1e-12)
    top <- which(one)[which.max(s$prob[one])]
    expect_equal(b$mode_cell[i], s$cell[top])
    expect_equal(b$mode_dist_m[i], d_mid[top], tolerance = 1e-9)
  }
  # Printed lines wrap at the console width; compare with the spaces collapsed.
  out <- gsub("\\s+", " ", paste(utils::capture.output(print(s)), collapse = " "))
  expect_match(out, "Misfit: 0 of 2 bins")
  # The radius is labelled as what it is, the farthest cell centre at p >= 0.5,
  # which falls short of the curve's d50 by up to a cell.
  expect_match(out, paste0("Near: 1 of 2 bins put under 0.001 of their probability within ",
                           "300 m of the receiver with the most detections \\(per receiver, ",
                           "its farthest cell centre at p >= 0.5\\)"))
  expect_match(out, "Modes sit a median [0-9,]+ m \\(max [0-9,]+ m\\) .* over 2 bins")
})

test_that("control$near_m sets the radius, and subsetting keeps the rule", {
  lake <- near_lake()
  ring <- counts_for(lake$field, c(1, 0, 0), 60, key = "station_id")
  s <- location_surface(ring, lake$field, "none", control = list(near_m = 1500, kept_mass = 1))
  b <- attr(s, "bins")
  expect_equal(b$near_m, 1500)
  xy <- terra::xyFromCell(lake$grid, s$cell)
  inside <- (xy[, 1] - 1525)^2 + (xy[, 2] - 1525)^2 <= 1500^2
  expect_equal(b$near_mass, sum(s$prob[inside]), tolerance = 1e-12)
  printed <- gsub("\\s+", " ", paste(utils::capture.output(print(s[s$prob > 0, ])),
                                     collapse = " "))
  expect_match(printed,
               "within 1500 m of the receiver with the most detections \\(control\\$near_m\\)")
  expect_error(location_surface(ring, lake$field, "none", control = list(near_m = -1)),
               "control\\$near_m")
  expect_error(location_surface(ring, lake$field, "none", control = list(near_m = Inf)),
               "control\\$near_m")
})

test_that("without a distance table the radius is NA and print says how to set it", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  s <- location_surface(counts_for(field, c(3, 1, 0, 0), 12), field, "none")
  b <- attr(s, "bins")
  expect_true(is.na(b$near_m) && is.na(b$near_mass))
  expect_equal(b$near_receiver, "R1")
  # The mode distance needs only positions, which the fixture holds.
  expect_false(is.na(b$mode_dist_m))
  expect_output(print(s), "Near: not computed; .* set control\\$near_m")
  given <- location_surface(counts_for(field, c(3, 1, 0, 0), 12), field, "none",
                            control = list(near_m = 60))
  expect_equal(attr(given, "bins")$near_m, 60)
  silent <- attr(location_surface(counts_for(field, c(0, 0, 0, 0), 12), field), "bins")
  expect_true(is.na(silent$near_receiver) && is.na(silent$mode_dist_m) &&
                is.na(silent$near_mass))
  bare <- field
  bare$receiver_xy <- NULL
  nb <- attr(location_surface(counts_for(field, c(3, 1, 0, 0), 12), bare, "none"), "bins")
  expect_true(is.na(nb$mode_dist_m) && is.na(nb$near_mass))
  expect_equal(nb$near_receiver, "R1")
})
