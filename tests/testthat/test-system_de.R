brute_at_least <- function(p, k) {
  s <- ncol(p)
  hits <- as.matrix(expand.grid(rep(list(0:1), s)))
  out <- numeric(nrow(p))
  for (h in seq_len(nrow(hits))) {
    if (sum(hits[h, ]) < k) next
    prob <- rep(1, nrow(p))
    for (j in seq_len(s)) prob <- prob * if (hits[h, j]) p[, j] else 1 - p[, j]
    out <- out + prob
  }
  out
}

test_that("at_least_k() matches enumeration over every detection pattern", {
  set.seed(1)
  p <- matrix(stats::runif(40), 10, 4)
  for (k in 1:4) expect_equal(at_least_k(p, k), brute_at_least(p, k))
  expect_equal(at_least_k(p, 5), rep(0, 10))
})

test_that("coverage per day uses the receivers listening that day", {
  grid <- peninsula_lake()
  d <- water_distances(grid, peninsula_receivers())
  f <- de_field(d, de_logistic(250, 600))
  listening <- data.frame(station_id = c("west", "east", "west", "east"),
                          date = as.Date(c("2024-06-01", "2024-06-01", "2024-06-02",
                                           "2024-06-03")),
                          listening = c(TRUE, TRUE, TRUE, FALSE))
  cover <- system_de(f, listening)
  expect_identical(names(cover), c("date", "cell", "p", "n_listening"))
  expect_equal(sort(unique(cover$date)), as.Date(c("2024-06-01", "2024-06-02")))
  both <- cover$p[cover$date == as.Date("2024-06-01")]
  expect_equal(both, 1 - (1 - f$p[, "west"]) * (1 - f$p[, "east"]))
  expect_equal(cover$p[cover$date == as.Date("2024-06-02")], unname(f$p[, "west"]))
  expect_equal(unique(cover$n_listening), c(2L, 1L))
  two <- system_de(f, listening, k = 2)
  expect_equal(two$p[two$date == as.Date("2024-06-01")], unname(f$p[, "west"] * f$p[, "east"]))
  r <- system_de(f, listening, format = "raster")
  expect_s4_class(r, "SpatRaster")
  expect_identical(names(r), c("2024-06-01", "2024-06-02"))
  expect_equal(terra::values(r[[1]], mat = FALSE)[f$cells], both)
})

test_that("serials map to stations through the deployment lookup", {
  grid <- peninsula_lake()
  rx <- data.frame(station_id = c("west", "east", "east"), receiver_serial = c("A", "B", "A"),
                   x = c(285, 885, 885), y = c(435, 435, 435),
                   deploy_datetime = as.POSIXct(c("2024-06-01 10:00", "2024-06-01 10:00",
                                                  "2024-06-05 15:00"), tz = "UTC"),
                   recover_datetime = as.POSIXct(c("2024-06-05 14:00", "2024-06-05 14:00", NA),
                                                 tz = "UTC"))
  f <- de_field(water_distances(grid, rx), de_logistic(250, 600))
  # On the swap day, 5 June, serial A belongs to its later deployment at east.
  listening <- data.frame(receiver_serial = c("A", "A", "A"),
                          date = as.Date(c("2024-06-04", "2024-06-05", "2024-06-30")))
  cover <- system_de(f, listening)
  expect_equal(cover$p[cover$date == as.Date("2024-06-04")], unname(f$p[, "west"]))
  expect_equal(cover$p[cover$date == as.Date("2024-06-05")], unname(f$p[, "east"]))
  expect_equal(cover$p[cover$date == as.Date("2024-06-30")], unname(f$p[, "east"]))
  undated <- rx[1:2, c("station_id", "receiver_serial", "x", "y")]
  g <- de_field(water_distances(grid, undated), de_logistic(250, 600))
  cover2 <- system_de(g, data.frame(receiver_serial = "B", date = as.Date("2024-06-04")))
  expect_equal(cover2$p, unname(g$p[, "east"]))
  moved <- rx[, c("station_id", "receiver_serial", "x", "y")]
  h <- de_field(water_distances(grid, moved), de_logistic(250, 600))
  expect_error(system_de(h, data.frame(receiver_serial = "A", date = as.Date("2024-06-04"))),
               "more than one station")
})

test_that("bad listening tables stop", {
  f <- de_field(water_distances(peninsula_lake(), peninsula_receivers()), de_logistic(250, 600))
  day <- as.Date("2024-06-01")
  expect_error(system_de(f, data.frame(station_id = "west", date = "2024-06-01")), "class Date")
  expect_error(system_de(f, data.frame(station_id = "north", date = day)),
               "names a receiver the field does not hold: north")
  expect_error(system_de(f, data.frame(serial = "west", date = day)), "needs a `station_id`")
  expect_error(system_de(f, data.frame(station_id = "west", date = day, listening = "yes")),
               "must be logical")
  expect_error(system_de(f, data.frame(station_id = "west", date = day, listening = FALSE)),
               "no listening rows")
  expect_error(system_de(f, data.frame(station_id = "west", date = day), k = 0), "whole number")
  expect_error(system_de(list(), data.frame(station_id = "west", date = day)), "de_field")
})

test_that("a final recovery day belongs to the deployment that ended", {
  # A listening row on a final recovery day has no later deployment of its
  # serial, so it goes to the station the serial left.
  grid <- peninsula_lake()
  rx <- data.frame(station_id = c("west", "east"), receiver_serial = c("A", "B"),
                   x = c(285, 885), y = c(435, 435),
                   deploy_datetime = as.POSIXct(c("2024-06-01 10:00", "2024-06-01 10:00"),
                                                tz = "UTC"),
                   recover_datetime = as.POSIXct(c("2024-06-05 14:00", NA), tz = "UTC"))
  f <- de_field(water_distances(grid, rx), de_logistic(250, 600))
  listening <- data.frame(receiver_serial = c("A", "A"),
                          date = as.Date(c("2024-06-04", "2024-06-05")))
  cover <- system_de(f, listening)
  expect_equal(cover$p[cover$date == as.Date("2024-06-05")], unname(f$p[, "west"]))
  by_station <- system_de(f, data.frame(station_id = "west", date = listening$date))
  expect_identical(cover, by_station)
  late <- data.frame(receiver_serial = c("A", "A", "B"),
                     date = as.Date(c("2024-06-06", "2024-06-07", "2024-05-01")))
  expect_error(system_de(f, late), paste("3 rows of `listening` name a receiver_serial",
                                         "with no deployment on that date: A on 2024-06-06"))
})

# ---- Values in 0 to 1 ----
test_that("every value lies in 0 to 1 and k past the receivers stops", {
  # Unclamped, 1 - rowSums(q) gives -4.4e-16 where no receiver hears a cell
  # (on a large grid at k = 3).
  p <- matrix(c(1e-300, 0.3, 1e-300, 1e-17, 0.2, 1e-300), 3, 2)
  for (k in 1:2) expect_true(all(at_least_k(p, k) >= 0 & at_least_k(p, k) <= 1))
  grid <- peninsula_lake()
  f <- de_field(water_distances(grid, peninsula_receivers()), de_logistic(250, 600))
  listening <- data.frame(station_id = c("west", "east", "west"),
                          date = as.Date(c("2024-06-01", "2024-06-01", "2024-06-02")))
  for (k in 1:2) {
    cover <- system_de(f, listening, k = k)
    expect_true(all(cover$p >= 0 & cover$p <= 1))
  }
  two <- system_de(f, listening, k = 2)
  expect_true(all(two$p[two$date == as.Date("2024-06-02")] == 0))
  expect_error(system_de(f, listening, k = 3), "the field holds 2 receivers")
})

# ---- Listening rows and dates ----
test_that("a field scored against listening rows of another date warns", {
  grid <- peninsula_lake()
  model <- de_model(de_logistic(250, 600), modulators = list(wind = function(date) 1))
  f <- de_field(water_distances(grid, peninsula_receivers()), model,
                date = as.Date("2024-06-02"))
  wrong_day <- data.frame(station_id = c("west", "east"), date = as.Date("2024-06-01"))
  expect_warning(cover <- system_de(f, wrong_day), "built for 2024-06-02")
  expect_identical(unique(cover$date), as.Date("2024-06-01"))
  right_day <- data.frame(station_id = c("west", "east"), date = as.Date("2024-06-02"))
  expect_no_warning(system_de(f, right_day))
  undated <- de_field(water_distances(grid, peninsula_receivers()), de_logistic(250, 600))
  expect_no_warning(system_de(undated, wrong_day))
})

test_that("the listening list of detqc's hand-off is read through its days table", {
  grid <- peninsula_lake()
  f <- de_field(water_distances(grid, peninsula_receivers()), de_logistic(250, 600))
  days <- data.frame(station_id = c("west", "east"), date = as.Date("2024-06-01"))
  handoff <- list(deployments = data.frame(station_id = c("west", "east")), days = days)
  expect_identical(system_de(f, handoff), system_de(f, days))
  expect_error(system_de(f, list(deployments = days)), "data frame with a `date` column")
})

# ---- Detection per bin ----
three_receiver_field <- function() {
  grid <- peninsula_lake()
  rx <- data.frame(station_id = c("west", "east", "mid"), x = c(285, 885, 585),
                   y = c(435, 435, 135))
  suppressMessages(de_field(water_distances(grid, rx), de_logistic(250, 600)))
}

test_that("transmissions = 1 is identical to the per-transmission result", {
  f <- three_receiver_field()
  listening <- data.frame(station_id = c("west", "east", "mid", "west"),
                          date = as.Date(c(rep("2024-06-01", 3), "2024-06-02")))
  for (k in 1:3) {
    expect_identical(system_de(f, listening, k = k, transmissions = 1),
                     system_de(f, listening, k = k))
  }
  expect_identical(system_de(f, listening, transmissions = 1, format = "raster")[],
                   system_de(f, listening, format = "raster")[])
  # The default output has the columns and no attributes beyond a data frame's.
  cover <- system_de(f, listening)
  expect_identical(names(attributes(cover)), c("names", "class", "row.names"))
})

test_that("with k = 1 a bin is heard with probability 1 - (1 - p)^n", {
  f <- three_receiver_field()
  listening <- data.frame(station_id = c("west", "east", "mid"),
                          date = as.Date("2024-06-01"))
  one <- system_de(f, listening)$p
  for (n in c(2, 24)) {
    expect_equal(system_de(f, listening, transmissions = n)$p, 1 - (1 - one)^n)
  }
  # One receiver alone: its own per-bin probability.
  west <- system_de(f, listening[1, ], transmissions = 24)$p
  expect_equal(west, 1 - (1 - unname(f$p[, "west"]))^24)
})

test_that("with k = 2 the per-bin value agrees with a Monte Carlo count", {
  f <- three_receiver_field()
  listening <- data.frame(station_id = c("west", "east", "mid"),
                          date = as.Date("2024-06-01"))
  n <- 6
  exact <- system_de(f, listening, k = 2, transmissions = n)$p
  # Simulate every transmission at every receiver for 20,000 bins per cell.
  cells <- unique(round(seq(1, length(f$cells), length.out = 12)))
  reps <- 20000
  set.seed(11)
  mc <- vapply(cells, function(i) {
    heard <- vapply(seq_len(3), function(s) {
      stats::rbinom(reps, n, f$p[i, s]) > 0
    }, logical(reps))
    mean(rowSums(heard) >= 2)
  }, numeric(1))
  se <- sqrt(pmax(exact[cells] * (1 - exact[cells]), 1 / reps) / reps)
  expect_true(all(abs(mc - exact[cells]) <= 3 * se),
              info = paste(signif(abs(mc - exact[cells]) / se, 3), collapse = " "))
  expect_gt(max(exact[cells]), 0.1)
  expect_lt(min(exact[cells]), 0.9)
})

test_that("a bad transmissions stops", {
  f <- three_receiver_field()
  day <- data.frame(station_id = "west", date = as.Date("2024-06-01"))
  for (bad in list(0, 1.5, c(2, 3), NA_real_, "24", Inf)) {
    expect_error(system_de(f, day, transmissions = bad),
                 "`transmissions` must be one whole number of at least 1.", fixed = TRUE)
  }
})

test_that("a day on which every receiver is off or deaf gives 0 in every cell", {
  f <- de_field(water_distances(peninsula_lake(), peninsula_receivers()), de_logistic(250, 600))
  days <- as.Date(c("2024-06-01", "2024-06-02", "2024-06-03"))
  listening <- data.frame(station_id = rep(c("west", "east"), 3), date = rep(days, each = 2),
                          listening = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE),
                          status = c("listening", "listening", "off", "deaf",
                                     "listening", "listening"))
  cover <- system_de(f, listening)
  expect_equal(sort(unique(cover$date)), days)
  off <- cover[cover$date == days[2], ]
  expect_equal(nrow(off), length(f$cells))
  expect_true(all(off$p == 0))
  expect_true(all(off$n_listening == 0L))
  both <- 1 - (1 - f$p[, "west"]) * (1 - f$p[, "east"])
  expect_equal(mean(cover$p), mean(c(both, rep(0, length(f$cells)), both)))
  expect_true(all(system_de(f, listening, transmissions = 24)$p[cover$date == days[2]] == 0))
  r <- system_de(f, listening, format = "raster")
  expect_identical(names(r), format(days))
  expect_true(all(terra::values(r[[2]], mat = FALSE)[f$cells] == 0))
  # Every day off gives zeros.
  all_off <- transform(listening, listening = FALSE, status = "off")
  expect_true(all(system_de(f, all_off)$p == 0))
})

test_that("a day whose receivers may have listened with their data missing is left out", {
  f <- de_field(water_distances(peninsula_lake(), peninsula_receivers()), de_logistic(250, 600))
  days <- as.Date(c("2024-06-01", "2024-06-02"))
  for (why in c("not_downloaded", "reference_silent", "silent_run")) {
    listening <- data.frame(station_id = rep(c("west", "east"), 2), date = rep(days, each = 2),
                            listening = c(TRUE, TRUE, FALSE, FALSE),
                            status = c("listening", "listening", "off", why))
    expect_equal(unique(system_de(f, listening)$date), days[1], info = why)
  }
  # Without a status column the reason is unknown, so the day is left out.
  plain <- data.frame(station_id = rep(c("west", "east"), 2), date = rep(days, each = 2),
                      listening = c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(unique(system_de(f, plain)$date), days[1])
  # NA listening is not a listening row and not an off one.
  unsure <- transform(plain, listening = c(TRUE, TRUE, NA, FALSE),
                      status = c("listening", "listening", "off", "off"))
  expect_equal(unique(system_de(f, unsure)$date), days[1])
  masked <- transform(plain, listening = FALSE, status = "not_downloaded")
  expect_error(system_de(f, masked), "no day on which every receiver is off or deaf")
})

test_that("on a mixed day a not_downloaded receiver is left out and the others count", {
  f <- de_field(water_distances(peninsula_lake(), peninsula_receivers()), de_logistic(250, 600))
  day <- data.frame(station_id = c("west", "east"), date = as.Date("2024-06-01"),
                    listening = c(TRUE, FALSE), status = c("listening", "not_downloaded"))
  cover <- system_de(f, day)
  expect_equal(cover$p, unname(f$p[, "west"]))
  expect_true(all(cover$n_listening == 1L))
})

test_that("a row naming a receiver the field lacks is an error on every kind of day", {
  f <- de_field(water_distances(peninsula_lake(), peninsula_receivers()), de_logistic(250, 600))
  days <- as.Date(c("2024-06-01", "2024-06-02"))
  listening <- data.frame(station_id = c("west", "east", "west", "north"),
                          date = rep(days, each = 2), listening = c(TRUE, TRUE, FALSE, FALSE),
                          status = c("listening", "listening", "off", "off"))
  expect_error(system_de(f, listening),
               "1 row of `listening` names a receiver the field does not hold: north.",
               fixed = TRUE)
  only_off <- listening[3:4, ]
  expect_error(system_de(f, only_off), "the field does not hold: north", fixed = TRUE)
  # A wrong receiver key is an error whatever the row's status.
  mixed <- data.frame(station_id = c("west", "north"), date = days[1],
                      listening = c(TRUE, FALSE), status = c("listening", "not_downloaded"))
  expect_error(system_de(f, mixed), "the field does not hold: north", fixed = TRUE)
  # A day left out of the coverage can still hold a wrong key.
  masked <- data.frame(station_id = c("west", "east", "north"), date = days[c(1, 2, 2)],
                       listening = c(TRUE, FALSE, FALSE),
                       status = c("listening", "not_downloaded", "not_downloaded"))
  expect_error(system_de(f, masked), "the field does not hold: north", fixed = TRUE)
  masked$station_id[3] <- "west"
  expect_equal(unique(system_de(f, masked)$date), days[1])
})
