# nolint start: object_usage_linter.
# The fixtures come from helper-fixtures.R, which the linter cannot see.

# ---- Known answers ----
test_that("the binomial surface equals brute-force Bayes on a 10 by 10 grid", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  for (case in list(list(y = c(5, 2, 0, 1), n = 12), list(y = c(0, 0, 0, 0), n = 30),
                    list(y = c(12, 12, 0, 3), n = 12))) {
    counts <- counts_for(field, case$y, case$n)
    got <- full_surface(counts, field, "none")
    want <- brute_posterior(field, case$y, case$n)
    expect_lt(max(abs(got - want$post)), 1e-10)
    s <- location_surface(counts, field, "none", control = list(kept_mass = 1))
    expect_lt(abs(attr(s, "bins")$log_lik - want$log_lik), 1e-10)
  }
})

test_that("quadrature agrees with Monte Carlo over 100,000 draws on one bin", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80,
                         sigma = 1.362)
  y <- c(6, 2, 0, 1)
  n <- 20
  got <- full_surface(counts_for(field, y, n), field)
  # Stratified Monte Carlo: one uniform draw per stratum of the normal, so
  # the 100,000 draws cover the distribution evenly.
  set.seed(42)
  z <- stats::qnorm((seq_len(1e5) - stats::runif(1e5)) / 1e5)
  ll <- numeric(nrow(field$eta))
  for (s in seq_along(y)) {
    for (c in seq_len(nrow(field$eta))) {
      ll[c] <- ll[c] + log(mean(stats::dbinom(y[s], n, stats::plogis(field$eta[c, s] + 1.362 * z))))
    }
  }
  mc <- exp(ll - max(ll))
  mc <- mc / sum(mc)
  expect_lt(max(abs(got - mc)), 1e-6)
})

test_that("the tabulated term matches a plain fixed-node logit-normal likelihood", {
  # The reference uses fixed nodes. Against a Riemann sum on this bin its
  # surface is 1.03e-06 off and the adaptive one 5.2e-15, so the bound here is
  # the reference's.
  field <- fixture_square_lake(sigma = 1)
  y <- c(9, 7, 3, 12, 3)
  got <- full_surface(counts_for(field, y, 20), field, nodes = 30)
  ll <- fixed_node_loglik_logitnormal(y, 20, field$eta, 1, nodes = 30)
  want <- exp(ll - max(ll))
  want <- want / sum(want)
  expect_lt(max(abs(got - want)), 2e-6)
})

test_that("the log-likelihood stays finite at 1,920 transmissions and no detection", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- counts_for(field, c(0, 0, 0, 0), 1920)
  for (form in c("station_hour", "none")) {
    s <- location_surface(counts, field, form, control = list(kept_mass = 1))
    expect_true(all(is.finite(s$prob)))
    expect_true(is.finite(attr(s, "bins")$log_lik))
    expect_equal(sum(s$prob), 1, tolerance = 1e-12)
  }
})

test_that("with sigma 0 the station-hour form is the binomial", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80,
                         sigma = 0)
  counts <- counts_for(field, c(5, 2, 0, 1), 12)
  expect_equal(full_surface(counts, field), full_surface(counts, field, "none"))
})

test_that("tables over the logit match direct adaptive quadrature", {
  od <- ls_overdispersion("station_hour", 1.362, ls_control(list()),
                          data.frame(bin_start = 0, bin_end = 3600))
  eta <- c(stats::runif(500, -40, 3), -Inf, Inf)
  prep <- ls_prepare(matrix(eta[is.finite(eta)], ncol = 1), od, 1)
  for (y in c(0, 1, 10, 30)) {
    got <- ls_term(eta, y, 30, prep)
    want <- ls_table(c(eta[1:500], 99), y, 30, od)$c0
    expect_lt(max(abs(got[1:500] - want)), 5e-9)
    expect_equal(got[501:502], c(if (y == 0) 0 else -Inf, if (y == 30) 0 else -Inf))
  }
})

test_that("adaptive nodes find the mode where fixed nodes miss it", {
  od <- list(sigma = 1.362, z = gauss_hermite(20)$z, w = gauss_hermite(20)$w)
  # Riemann sum on a fine grid of z around the mode as the reference.
  ref <- function(eta, y, n, s) {
    z <- ls_mode(eta, y, n, s) + seq(-40, 40, length.out = 160001)
    v <- lchoose(n, y) + y * stats::plogis(eta + s * z, log.p = TRUE) +
      (n - y) * stats::plogis(-(eta + s * z), log.p = TRUE) + stats::dnorm(z, log = TRUE)
    max(v) + log(sum(exp(v - max(v))) * (z[2] - z[1]))
  }
  for (eta in c(-8, -1, 2)) {
    got <- ls_table(c(eta, 99), 100, 200, od)$c0[1]
    expect_lt(abs(got - ref(eta, 100, 200, 1.362)), 1e-7)
  }
  expect_equal(gauss_hermite(1), list(z = 0, w = 1))
  expect_equal(sum(gauss_hermite(20)$w), 1)
})

# ---- Square-lake known-truth example ----
test_that("binomial form: the square lake's inside case gives a 95% region of 183 cells", {
  # 5,489 cells is the 95% region of positionR 1.7.0 on the same grid,
  # stations, curve and counts; positionR 1.7.0 is not run here.
  v1_cells <- 5489
  field <- fixture_square_lake()
  prob <- full_surface(counts_for(field, c(9, 7, 3, 12, 3), 20), field, "none")
  truth <- match(column_major_to_terra(6055), field$cells)
  expect_equal(length(region95(prob)), 183)
  expect_equal(prob[truth], 0.0165137105191619, tolerance = 1e-9)
  expect_true(truth %in% region95(prob))
  expect_lt(length(region95(prob)), v1_cells / 25)
  edge <- full_surface(counts_for(field, c(2, 3, 0, 9, 1), 20), field, "none")
  expect_equal(length(region95(edge)), 483)
  expect_equal(edge[match(column_major_to_terra(8163), field$cells)], 0.00180075546819677,
               tolerance = 1e-9)
})

test_that("binomial form: four far silent stations barely move the surface (183 to 178)", {
  # Silent receivers are evidence and must stay in.
  field <- fixture_square_lake(stations = 9)
  prob <- full_surface(counts_for(field, c(9, 7, 3, 12, 3, 0, 0, 0, 0), 20), field, "none")
  expect_equal(length(region95(prob)), 178)
  expect_equal(prob[match(column_major_to_terra(6055), field$cells)], 0.0171337944542997,
               tolerance = 1e-9)
})

test_that("binomial form only: more transmissions shrink the region (183 to 27 cells)", {
  y200 <- c(87, 43, 34, 81, 19)
  field <- fixture_square_lake()
  prob <- full_surface(counts_for(field, y200, 200), field, "none")
  expect_equal(length(region95(prob)), 27)
})

# ---- Safeguards ----
test_that("a degenerate surface stops with the animal and bin named", {
  # A surface with one value stops here, before a later script reads it.
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0 * d - 1)
  counts <- counts_for(field, c(1, 0, 0, 0), 10, animal = "pike-7", start = 5)
  expect_error(location_surface(counts, field, "none"),
               "Animal pike-7, bin starting 2024-06-01 05:00:00 UTC: .*same likelihood")
  field$eta[, 1] <- -Inf
  expect_error(location_surface(counts, field), "zero likelihood")
  expect_match(ls_degenerate(c(NaN, 1)), "NaN")
  expect_error(ls_stop_degenerate(ls_degenerate(c(NaN, 1)), "a",
                                  as.POSIXct("2024-01-01", tz = "UTC")),
               "NaN", class = "detfield_error_degenerate_bin")
  field$eta[] <- -Inf
  expect_error(location_surface(counts_for(field, c(0, 0, 0, 0), 10), field), "same likelihood")
  expect_null(ls_degenerate(5))
})

test_that("each bin keeps the smallest set of cells holding the kept mass", {
  # The output filter is the kept mass, stated.
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- counts_for(field, c(3, 1, 0, 0), 12)
  for (mass in c(0.5, 0.95, 0.999)) {
    s <- location_surface(counts, field, control = list(kept_mass = mass))
    expect_gte(sum(s$prob), mass * (1 - 1e-12))
    expect_lt(sum(s$prob) - min(s$prob), mass)
    expect_equal(attr(s, "bins")$kept_mass, sum(s$prob))
    expect_equal(attr(s, "kept")$mass, mass)
  }
})

test_that("bins are independent and one call holds many animals", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  a <- counts_for(field, c(3, 1, 0, 0), 12, "a", 0)
  b <- counts_for(field, c(0, 0, 2, 5), 12, "b", 3)
  c2 <- counts_for(field, c(0, 1, 1, 0), 30, "a", 1)
  both <- location_surface(rbind(a, b, c2), field)
  alone <- location_surface(b, field)
  expect_equal(both$prob[both$animal_id == "b"], alone$prob)
  expect_equal(nrow(attr(both, "bins")), 3)
  expect_equal(attr(both, "bins")$n_transmissions, c(12, 30, 12))
})

test_that("two cores give the same surfaces as one", {
  # Parallel runs go through control$cores.
  skip_on_os("windows")
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- rbind(counts_for(field, c(3, 1, 0, 0), 12, "a"),
                  counts_for(field, c(0, 0, 2, 5), 12, "b"))
  expect_equal(location_surface(counts, field, control = list(cores = 2)),
               location_surface(counts, field))
  field$eta[, 1] <- -Inf
  expect_error(location_surface(counts, field, control = list(cores = 2)), "Animal a")
})

test_that("a keyed data.table of counts gives the same surface", {
  # A stale data.table key can return partial joins silently.
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- counts_for(field, c(3, 1, 0, 0), 12)
  keyed <- data.table::as.data.table(counts[4:1, ])
  data.table::setkey(keyed, receiver_serial)
  expect_equal(location_surface(keyed, field)$prob, location_surface(counts, field)$prob)
})

test_that("the field's logits are used and its model is never called", {
  # Detection probability is computed once, in the field.
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  expect_error(field$model$predict(), "must not call")
  expect_s3_class(location_surface(counts_for(field, c(3, 1, 0, 0), 12), field), "loc_surface")
  p_only <- field
  p_only$eta <- NULL
  p_only$p <- stats::plogis(field$eta)
  expect_equal(full_surface(counts_for(field, c(3, 1, 0, 0), 12), p_only),
               full_surface(counts_for(field, c(3, 1, 0, 0), 12), field), tolerance = 1e-10)
})

test_that("receivers match the field by serial or by station", {
  rx <- fixture_receivers()
  field <- fixture_field(fixture_grid(10), rx, function(d) 0.5 - d / 80, key = "station_id")
  counts <- counts_for(field, c(3, 1, 0, 0), 12, key = "station_id")
  counts$receiver_serial <- paste0("VR2W-", 1:4)
  s <- location_surface(counts, field)
  expect_equal(attr(s, "key"), "station_id")
  field$key <- NULL
  expect_equal(location_surface(counts, field)$prob, s$prob)
  counts$station_id[2] <- "S9"
  expect_error(location_surface(counts, field),
               "1 listening receiver in `counts` has .* \\(matched by station_id\\): S9")
  # A key named by the field is checked too. A station with no column in
  # the field would otherwise enter the surface as NA.
  field$key <- "station_id"
  expect_error(location_surface(counts, field), "\\(matched by station_id\\): S9")
  counts$station_id <- NULL
  expect_error(location_surface(counts, field), "a column `counts` lacks")
})

test_that("rows of a bin must agree and must not repeat a receiver", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- counts_for(field, c(3, 1, 0, 0), 12)
  bad <- counts
  bad$n[2] <- 13L
  expect_error(location_surface(bad, field), "rows give different `n`")
  expect_error(location_surface(rbind(counts, counts[1, ]), field), "R1 appears twice")
  bad <- counts
  bad$y[1] <- 30L
  expect_message(s <- location_surface(bad, field), "raised .* in 1 of 1 bins")
  expect_equal(attr(s, "bins")$n_transmissions, 30)
  bad <- counts
  bad$y[1] <- NA
  expect_error(location_surface(bad, field), "whole numbers")
  bad <- counts
  bad$bin_end <- bad$bin_start
  expect_error(location_surface(bad, field), "must follow")
})

test_that("rows with listening FALSE are left out with a message", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- counts_for(field, c(3, 1, 0, 0), 12)
  counts$listening[4] <- FALSE
  msg <- paste(testthat::capture_messages(s <- location_surface(counts, field)), collapse = "")
  expect_match(msg, "1 row \\(0 detections\\) with listening = FALSE was left out\\.")
  # With no bin lost, the message names none.
  expect_no_match(msg, "with them")
  expect_equal(attr(s, "bins")$n_listening, 3)
  counts$listening <- FALSE
  expect_error(suppressMessages(location_surface(counts, field)), "no listening receiver")
})

test_that("overdispersion settings are checked", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- counts_for(field, c(3, 1, 0, 0), 12)
  day <- counts
  day$bin_end <- day$bin_start + 86400
  expect_warning(location_surface(day, field), "describes hour bins")
  expect_silent(location_surface(day, field, "none"))
  no_sigma <- field
  no_sigma$model$sigma$station_hour <- NA
  expect_error(location_surface(counts, no_sigma), "gives no `sigma\\$station_hour`")
  # The message gives the call that fixes it, on one line that can be copied
  expect_error(location_surface(counts, no_sigma),
               "\n  de_field(distances, de_model(model, sigma = list(station_hour = 0)))\n",
               fixed = TRUE)
  expect_s3_class(location_surface(counts, no_sigma, "none"), "loc_surface")
  s <- location_surface(counts, field)
  expect_equal(attr(s, "overdispersion"), list(form = "station_hour", sigma = 2.55, nodes = 64L))
  expect_warning(s <- location_surface(counts, field, control = list(nodes = 20)),
                 "below the 64 nodes")
  expect_equal(attr(s, "overdispersion")$nodes, 20L)
})

test_that("control and field are checked", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- counts_for(field, c(3, 1, 0, 0), 12)
  expect_error(location_surface(counts, field, control = list(node = 5)), "Unknown setting")
  expect_error(location_surface(counts, field, control = list(kept_mass = 2)), "kept_mass")
  expect_error(location_surface(counts, field, control = list(cores = 0)), "cores")
  expect_error(location_surface(counts, field, control = "x"), "must be a list")
  expect_error(location_surface(counts, list(eta = "x")), "numeric matrix")
  expect_error(location_surface(counts, "x"), "must be a detection field")
  bad <- field
  bad$cells <- bad$cells[-1]
  expect_error(location_surface(counts, bad), "one row per cell")
})

test_that("the result carries what a reader needs to check it", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  s <- location_surface(counts_for(field, c(3, 1, 0, 0), 12), field)
  expect_named(s, c("animal_id", "bin_start", "cell", "prob"))
  b <- attr(s, "bins")
  expect_named(b, c("animal_id", "bin_start", "bin_end", "n_transmissions", "n_listening",
                    "n_detecting", "y_total", "log_lik", "kept_mass", "n_cells", "misfit_p",
                    "misfit_receiver", "mode_cell", "mode_dist_m", "near_receiver", "near_m",
                    "near_mass"))
  expect_equal(c(b$n_listening, b$n_detecting, b$y_total), c(4, 2, 4))
  expect_equal(attr(s, "n_water_cells"), 100L)
  expect_s4_class(sl_grid(s), "SpatRaster")
  expect_equal(nrow(attr(s, "receiver_xy")), 4)
  bare <- field
  bare$receiver_xy <- NULL
  bare$grid <- NULL
  s_bare <- location_surface(counts_for(field, c(3, 1, 0, 0), 12), bare)
  expect_null(attr(s_bare, "receiver_xy"))
  expect_null(sl_grid(s_bare))
  bare$grid <- "not/a/file.tif"
  expect_error(location_surface(counts_for(field, c(3, 1, 0, 0), 12), bare),
               "grid file not/a/file.tif does not exist")
  bare$grid <- 5
  expect_error(location_surface(counts_for(field, c(3, 1, 0, 0), 12), bare),
               "must be a SpatRaster")
  expect_output(print(s), "1 animal, 1 bin, 100 cells kept in all")
  tmp <- tempfile(fileext = ".rds")
  saveRDS(s, tmp)
  expect_s4_class(sl_grid(readRDS(tmp)), "SpatRaster")
})

# ---- The station-hour form, numerics and messages ----
test_that("under the station-hour form the square lake's regions are far larger than 183 cells", {
  # 183, 178 and 27 hold for the binomial only. These sizes match an
  # independent posterior built on integrate()
  # in all 16 case-sigma pairs.
  cells <- function(stations, y, n, sigma) {
    field <- fixture_square_lake(stations = stations, sigma = sigma)
    length(region95(full_surface(counts_for(field, y, n), field)))
  }
  inside <- c(9, 7, 3, 12, 3)
  expect_equal(cells(5, inside, 20, 1.362), 1582)
  expect_equal(cells(5, inside, 20, 2.55), 4642)
  expect_equal(cells(5, inside, 20, 2.66), 4939)
  # At sigma 2.55 the edge case is larger than positionR 1.7.0's 5,383 cells,
  # and ten times the transmissions leave it at 5,168 cells.
  expect_equal(cells(5, c(2, 3, 0, 9, 1), 20, 2.55), 7210)
  expect_equal(cells(5, c(87, 43, 34, 81, 19), 200, 2.55), 5168)
  expect_equal(cells(9, c(inside, 0, 0, 0, 0), 20, 2.55), 4342)
})

test_that("the default nodes grow with sigma and hold the term near its exact value", {
  # 20 nodes were 1.8e-2 logits off one receiver's term at sigma 6.
  expect_equal(vapply(c(0.5, 1.362, 2.55, 3.4, 6, 10), ls_nodes, integer(1)),
               c(20L, 35L, 64L, 85L, 150L, 200L))
  exact <- function(e, y, n, s) {
    h <- function(z) {
      lchoose(n, y) + y * stats::plogis(e + s * z, log.p = TRUE) +
        (n - y) * stats::plogis(-(e + s * z), log.p = TRUE) + stats::dnorm(z, log = TRUE)
    }
    zm <- stats::optimize(h, c(-60, 60), maximum = TRUE)$maximum
    hm <- h(zm)
    hm + log(stats::integrate(function(z) exp(h(z) - hm), zm - 30, zm + 30,
                              rel.tol = 1e-12, subdivisions = 2000)$value)
  }
  for (s in c(2.55, 3.4)) {
    gh <- gauss_hermite(ls_nodes(s))
    od <- list(sigma = s, z = gh$z, w = gh$w)
    for (y in c(0, 1, 15, 30)) {
      e <- c(-12, -6, -2, 0, 2, 99)
      got <- ls_table(e, y, 30, od)$c0
      want <- vapply(e[-6], exact, numeric(1), y = y, n = 30, s = s)
      expect_lt(max(abs(got - want)), 1e-5)
    }
  }
})

test_that("a very negative logit costs no more than an ordinary one", {
  # Without the closed form one cell at -10,000 logits costs seconds and
  # gigabytes, because the term table spans every logit of the field.
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- counts_for(field, c(3, 1, 0, 0), 12)
  base <- full_surface(counts, field)
  field$eta[1, ] <- -50000
  t <- system.time(far <- full_surface(counts, field))[["elapsed"]]
  expect_lt(t, 2)
  expect_equal(far[1], 0)
  expect_lt(max(abs(far[-1] / sum(far[-1]) - base[-1] / sum(base[-1]))), 1e-6)
  # The closed form below the table meets the table at its lower end.
  od <- ls_overdispersion("station_hour", 2.55, ls_control(list()),
                          data.frame(bin_start = 0, bin_end = 3600))
  prep <- ls_prepare(matrix(c(-200, 3), ncol = 1), od, 1)
  for (y in c(0, 1, 5)) {
    lo <- ls_lo(prep, y)
    edge <- ls_table(c(lo, lo + 0.01), y, 30, od)$c0[1]
    expect_lt(abs(ls_low(lo, y, 30, 2.55) - edge), 1e-9 * max(1, abs(edge)))
    expect_equal(ls_term(c(lo - 1, -Inf), y, 30, prep),
                 c(ls_low(lo - 1, y, 30, 2.55), if (y == 0) 0 else -Inf))
  }
})

test_that("a degenerate bin can be skipped and recorded instead of stopping the run", {
  # Skipping keeps the other bins of the run.
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- rbind(counts_for(field, c(3, 1, 0, 0), 12, "a1"),
                  counts_for(field, c(0, 2, 0, 0), 12, "a2"),
                  counts_for(field, c(4, 0, 0, 0), 12, "a3", start = 2))
  field$eta[, 1] <- -Inf
  expect_error(location_surface(counts, field), "Animal a1, .*on_degenerate")
  for (cores in c(1, 2)) {
    expect_warning(s <- location_surface(counts, field,
                                         control = list(on_degenerate = "skip", cores = cores)),
                   "2 of 3 bins were degenerate")
    expect_equal(unique(s$animal_id), "a2")
    skipped <- attr(s, "skipped")
    expect_equal(skipped$animal_id, c("a1", "a3"))
    expect_match(skipped$reason, "zero likelihood")
    expect_equal(nrow(attr(s, "bins")), 1)
  }
  expect_output(print(s), "Left out: 2 degenerate bins")
  expect_error(location_surface(counts[counts$animal_id != "a2", ], field,
                                control = list(on_degenerate = "skip")),
               "^All 2 bins are degenerate", class = "detfield_error_degenerate")
  # A stop on one bin has its own class, so a caller that catches the run with
  # every bin degenerate does not catch it; with two cores the class survives
  # the worker.
  for (cores in c(1, 2)) {
    err <- tryCatch(location_surface(counts, field, control = list(cores = cores)),
                    error = function(e) e)
    expect_s3_class(err, c("detfield_error_degenerate_bin", "detfield_error"))
    expect_false(inherits(err, "detfield_error_degenerate"))
  }
  expect_error(location_surface(counts, field, control = list(on_degenerate = "ignore")),
               "on_degenerate")
})

test_that("each bin reports how well the surface explains its detecting receivers", {
  # The second bin has 12 detections at each of two opposite receivers, so
  # the surface puts almost no mass near either and misfit_p flags it.
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  fits <- counts_for(field, c(10, 3, 3, 0), 12)
  clash <- counts_for(field, c(12, 0, 0, 12), 12, start = 1)
  s <- location_surface(rbind(fits, clash), field, "none", control = list(kept_mass = 1))
  b <- attr(s, "bins")
  expect_gt(b$misfit_p[1], 0.05)
  expect_lt(b$misfit_p[2], 1e-3)
  expect_true(b$misfit_receiver[2] %in% c("R1", "R4"))
  # Under the binomial it is the posterior predictive tail, by brute force.
  post <- brute_posterior(field, c(12, 0, 0, 12), 12)$post
  p <- stats::plogis(field$eta)
  tails <- vapply(c(1, 4), function(k) sum(post * (p[, k]^12)), numeric(1))
  expect_equal(b$misfit_p[2], min(tails), tolerance = 1e-10)
  # Under the station-hour form the tail is integrated over the effect.
  sh <- attr(location_surface(clash, field, control = list(kept_mass = 1)), "bins")
  post <- full_surface(clash, field)
  z <- stats::qnorm((seq_len(2e4) - 0.5) / 2e4)
  tail_sh <- vapply(c(1, 4), function(k) {
    sum(post * vapply(field$eta[, k], function(e) mean(stats::plogis(e + 2.55 * z)^12),
                      numeric(1)))
  }, numeric(1))
  expect_equal(sh$misfit_p, min(tail_sh), tolerance = 1e-3)
  expect_gt(sh$misfit_p, b$misfit_p[2])
  silent <- location_surface(counts_for(field, c(0, 0, 0, 0), 12), field)
  expect_true(is.na(attr(silent, "bins")$misfit_p))
  expect_output(print(s), "Misfit: 1 of 2 bins")
})

test_that("a grid given as a relative path is stored with the surface", {
  # A relative path read from another working directory would give NA
  # coordinates and areas without a warning.
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  dir <- tempfile("grid")
  dir.create(dir)
  terra::writeRaster(field$grid, file.path(dir, "g.tif"))
  field$grid <- "g.tif"
  s <- withr::with_dir(dir, location_surface(counts_for(field, c(3, 1, 0, 0), 12), field))
  expect_s4_class(attr(s, "grid"), "PackedSpatRaster")
  expect_false(anyNA(sample_locations(s, 5, seed = 1)$x))
  expect_false(anyNA(attr(credible_region(s), "summary")$area_m2))
  attr(s, "grid") <- "g.tif"
  expect_error(sample_locations(s, 5), "grid \\(g.tif\\) cannot be read")
  expect_warning(r <- credible_region(s), "cannot be read")
  expect_true(all(is.na(attr(r, "summary")$area_m2)))
})

test_that("field cells on land, mixed bin ends and n = 0 stop with the cause named", {
  # Each of the three stops with its cause.
  rx <- fixture_receivers()
  field <- fixture_field(fixture_grid(10), rx, function(d) 0.5 - d / 80)
  counts <- counts_for(field, c(3, 1, 0, 0), 12)
  land <- field
  land$grid <- fixture_grid(10, land = 1:4)
  expect_error(location_surface(counts, land), "4 of the field's 100 cells are not water")
  # A template grid without values cannot be checked and is accepted.
  template <- field
  template$grid <- terra::rast(field$grid)
  expect_s3_class(location_surface(counts, template), "loc_surface")
  template$grid <- suppressWarnings(terra::wrap(template$grid))
  expect_s3_class(location_surface(counts, template), "loc_surface")
  mixed <- counts
  mixed$bin_end[2] <- mixed$bin_end[2] + 1800
  expect_error(location_surface(mixed, field, "none"), "rows give different `bin_end`")
  zero <- counts_for(field, c(0, 0, 0, 0), 0)
  expect_error(location_surface(zero, field), "n is 0")
})

test_that("a bin whose rows all have listening FALSE is counted as lost", {
  # Messages wrap at the console width; a wide console keeps this one on a line.
  withr::local_options(width = 200)
  # The message counts the bins left out as well as the rows.
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  gone <- counts_for(field, c(3, 1, 0, 0), 12, start = 1)
  gone$listening <- FALSE
  expect_message(s <- location_surface(rbind(counts_for(field, c(3, 1, 0, 0), 12), gone), field),
                 "4 rows \\(4 detections\\) .* 1 bin with no listening receiver")
  expect_equal(nrow(attr(s, "bins")), 1)
})

test_that("a field whose own `$` method computes elements is read as stored", {
  # A `$` method that builds `p` from `eta` must not be called, so a field
  # holding only `p` must be read as stored.
  registerS3method("$", "trap_field", function(x, name) stop("`$` was called"))
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- counts_for(field, c(3, 1, 0, 0), 12)
  want <- location_surface(counts, field, control = list(kept_mass = 1))
  p_only <- unclass(field)
  p_only$p <- stats::plogis(p_only$eta)
  p_only$eta <- NULL
  class(p_only) <- c("trap_field", "de_field")
  got <- location_surface(counts, p_only, control = list(kept_mass = 1))
  expect_equal(got$prob[order(got$cell)], want$prob[order(want$cell)], tolerance = 1e-10)
})

# nolint end

# ---- Help pages, node counts and large runs ----
rd_text <- function(topic) {
  src <- system.file("man", paste0(topic, ".Rd"), package = "detfield")
  rd <- if (nzchar(src)) tools::parse_Rd(src) else tools::Rd_db("detfield")[[paste0(topic, ".Rd")]]
  gsub("\\s+", " ", paste(utils::capture.output(tools::Rd2txt(rd)), collapse = " "))
}

test_that("a station-day spread stops instead of being ignored", {
  m <- de_logistic(250, 600)
  expect_error(de_model(m, sigma = list(station_hour = 1, station_day = 3)),
               "station_day` is not used")
  expect_silent(de_model(m, sigma = list(station_hour = 1, station_day = NA)))
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  field$model$sigma$station_day <- 3
  expect_error(location_surface(counts_for(field, c(3, 1, 0, 0), 12), field),
               "station_day` is not used")
})

test_that("sigma is tied to its curve, with no bare value recommended", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  field$model$sigma$station_hour <- NA
  msg <- tryCatch(location_surface(counts_for(field, c(3, 1, 0, 0), 12), field),
                  error = conditionMessage)
  # The message is wrapped to the console width, so spaces are compared as one
  expect_match(gsub("\\s+", " ", msg), "estimated around the same curve")
  for (topic in c("location_surface", "de_model")) {
    txt <- rd_text(topic)
    expect_match(txt, "curve it was estimated", fixed = TRUE)
  }
})

test_that("too few nodes and a spread past the node cap warn", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- counts_for(field, c(3, 1, 0, 0), 12)
  expect_warning(location_surface(counts, field, control = list(nodes = 1)), "below the 64")
  field$model$sigma$station_hour <- 9
  expect_warning(location_surface(counts, field), "200 is the most allowed")
})

test_that("a spread of 0 is labelled as the binomial", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- counts_for(field, c(3, 1, 0, 0), 12)
  field$model$sigma$station_hour <- 0
  s0 <- location_surface(counts, field)
  s_none <- location_surface(counts, field, "none")
  expect_equal(attr(s0, "overdispersion"), list(form = "none", sigma = 0, nodes = NA_integer_))
  expect_equal(s0$prob, s_none$prob)
})

test_that("a bin heard only at receivers not listening is left out of the surface", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 2 - d / 60)
  heard <- counts_for(field, c(9, 4, 0, 0), 30, start = 0)
  lost <- counts_for(field, c(9, 4, 0, 0), 30, start = 1)
  lost$listening[1:2] <- FALSE
  counts <- rbind(heard, lost)
  expect_warning(s <- suppressMessages(location_surface(counts, field)),
                 "1 bin had detections only at receivers not listening \\(13 detections\\)")
  expect_equal(nrow(attr(s, "bins")), 1L)
  expect_equal(attr(s, "not_listening_bins")$detections, 13L)
  expect_output(print(s), "1 bin heard only at receivers not listening")
  expect_warning(suppressMessages(location_surface(counts, field)),
                 "without them the bin would read as unheard, so it is left out")
  lost2 <- counts_for(field, c(9, 4, 0, 0), 30, start = 2)
  lost2$listening[1:2] <- FALSE
  expect_warning(suppressMessages(location_surface(rbind(counts, lost2), field)),
                 paste("2 bins had detections only at receivers not listening \\(26 detections\\);",
                       "without them the bins would read as unheard, so they are left out"))
  none <- counts
  none$listening <- FALSE
  expect_error(suppressMessages(location_surface(none, field)),
               "all 8 rows have listening = FALSE. Check that the listening table")
})

test_that("a run past one table's rows stops with the remedy named", {
  expect_error(ls_check_rows(3221225472), "Set control\\$write_dir")
  expect_silent(ls_check_rows(2e9))
  expect_error(ls_check_rows(11, limit = 10), "11 rows, more than one table can hold \\(10\\)")
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- counts_for(field, c(3, 1, 0, 0), 12)
  local_mocked_bindings(ls_max_rows = function() 50)
  expect_message(s <- location_surface(counts, field, control = list(kept_mass = 0.5)),
                 "could keep more rows than one table holds")
  expect_lt(nrow(s), 50)
  expect_error(suppressMessages(location_surface(counts, field, control = list(kept_mass = 1))),
               "100 rows, more than one table can hold \\(50\\)")
})

test_that("write_dir saves each animal's surfaces and returns a table of files", {
  field <- fixture_field(fixture_grid(10), fixture_receivers(), function(d) 0.5 - d / 80)
  counts <- rbind(counts_for(field, c(3, 1, 0, 0), 12, "fish/1"),
                  counts_for(field, c(0, 2, 5, 0), 12, "fish 2", start = 1),
                  counts_for(field, c(1, 0, 0, 4), 12, "fish 2", start = 2))
  whole <- location_surface(counts, field)
  dir <- withr::local_tempdir()
  for (cores in c(1, 2)) {
    files <- location_surface(counts, field, control = list(write_dir = dir, cores = cores))
    expect_setequal(files$animal_id, c("fish 2", "fish/1"))
    expect_equal(files$bins[match(c("fish/1", "fish 2"), files$animal_id)], c(1L, 2L))
    expect_true(all(file.exists(files$file)))
    expect_equal(nrow(attr(files, "bins")), 3L)
    back <- lapply(files$file, readRDS)
    expect_s3_class(back[[1]], "loc_surface")
    expect_equal(sum(vapply(back, nrow, integer(1))), nrow(whole))
    one <- back[[match("fish/1", files$animal_id)]]
    wb <- attr(whole, "bins")
    expect_equal(one$prob, whole$prob[whole$animal_id == "fish/1"])
    expect_equal(attr(one, "bins")$log_lik, wb$log_lik[wb$animal_id == "fish/1"])
  }
  expect_error(location_surface(counts, field, control = list(write_dir = 3)), "write_dir")
  blocker <- withr::local_tempfile()
  writeLines("a file where the folder goes", blocker)
  inside_file <- list(write_dir = file.path(blocker, "x"))
  expect_error(suppressWarnings(location_surface(counts, field, control = inside_file)),
               "could not be created")
})

test_that("print gives the cells kept and their median per bin", {
  # Printed lines wrap at the console width; a wide console keeps each on a line.
  withr::local_options(width = 200)
  grid <- terra::rast(system.file("extdata", "stoney_depth.tif", package = "detfield"))
  field <- suppressMessages(de_field(water_distances(grid, stoney_track$receivers),
                                     de_model(de_logistic(c(250, 750), c(500, 1500),
                                                          depth = c(1, 35)),
                                              sigma = list(station_hour = 1))))
  one <- stoney_track$detections[stoney_track$detections$animal_id == "stoney-01", ]
  counts <- bin_detections(one, stoney_track$receivers, stoney_track$tags)
  s <- location_surface(counts, field)
  out <- paste(utils::capture.output(print(s)), collapse = "\n")
  n_cells <- attr(s, "bins")$n_cells
  # stoney-01 is detected in 43 hours.
  expect_match(out, sprintf("1 animal, 43 bins, %s cells kept in all",
                            format(nrow(s), big.mark = ",")))
  expect_match(out, sprintf("median %s of 4,818 water cells \\(range %d to %d\\)",
                            format(stats::median(n_cells)), min(n_cells), max(n_cells)))
  expect_match(out, "Spread: 1 logit per receiver-hour \\(sigma\\$station_hour\\)")
  expect_match(out, "integrated over [0-9]+ quadrature points")
  expect_false(grepl("cells kept of", out))
})

test_that("counts keyed by receiver_serial position on a field keyed by station", {
  grid <- terra::rast(system.file("extdata", "stoney_depth.tif", package = "detfield"))
  rx <- stoney_track$receivers
  field <- suppressMessages(de_field(water_distances(grid, rx),
                                     de_logistic(c(250, 750), c(500, 1500), depth = c(1, 35))))
  expect_identical(field$key, "station_id")
  det <- stoney_track$detections
  # Every receiver listens on both days of the study
  days <- expand.grid(receiver_serial = rx$receiver_serial,
                      date = as.Date(c("2024-06-01", "2024-06-02")), stringsAsFactors = FALSE)
  by_serial <- bin_detections(det, days, stoney_track$tags)
  expect_false("station_id" %in% names(by_serial))
  s1 <- location_surface(by_serial, field, overdispersion = "none")
  by_station <- by_serial
  by_station$station_id <- rx$station_id[match(by_station$receiver_serial, rx$receiver_serial)]
  s2 <- location_surface(by_station, field, overdispersion = "none")
  expect_equal(s1$prob, s2$prob)
  bad <- by_serial
  bad$receiver_serial[1] <- "999999"
  expect_error(location_surface(bad, field, overdispersion = "none"), "no deployment in the field")
})
