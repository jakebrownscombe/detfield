# nolint start: object_usage_linter.
# The fixtures come from helper-fixtures.R, which the linter cannot see.

# ---- Coverage of the 95% region on the square lake ----
# Heavy: 1,000 to 20,000 surfaces each. Skipped on CRAN; seeds fixed.

replay_draws <- function(field_eta, n_rep = 500) {
  # 500 detected bins per scenario from seed 7, scenario A first: a cell drawn
  # at random, five counts out of 20, and a draw with no detection skipped.
  set.seed(7)
  out <- list()
  for (scn in c("A", "B")) {
    kept <- 0
    while (kept < n_rep) {
      tc <- sample.int(10000, 1)
      y <- if (scn == "A") {
        stats::rbinom(5, 20, stats::plogis(field_eta[tc, ]))
      } else {
        stats::rbinom(5, 20, stats::plogis(field_eta[tc, ] + stats::rnorm(5, 0, 1)))
      }
      if (sum(y) == 0) next
      kept <- kept + 1
      out[[length(out) + 1]] <- list(scn = scn, tc = tc, y = y)
    }
  }
  out
}

coverage_of <- function(draws, field, form, nodes = 20) {
  counts <- do.call(rbind, lapply(seq_along(draws), function(i) {
    counts_for(field, draws[[i]]$y, 20, animal = "rep", start = i)
  }))
  s <- location_surface(counts, field, form, control = list(kept_mass = 0.95, nodes = nodes))
  key <- paste(s$bin_start, s$cell)
  start <- as.POSIXct("2024-06-01", tz = "UTC") + seq_along(draws) * 3600
  truth <- vapply(draws, function(d) column_major_to_terra(d$tc), numeric(1))
  list(coverage = mean(paste(start, truth) %in% key),
       median_cells = stats::median(attr(s, "bins")$n_cells))
}

test_that("500-bin coverage on the square lake holds under scenarios A and B", {
  skip_on_cran()
  field <- fixture_square_lake(sigma = 1)
  # The field lists cells in terra's order and the draws number cells column
  # by column; column_major_cells puts the logits in that order for the draws.
  column_major_cells <- match(column_major_to_terra(1:10000), field$cells)
  draws <- replay_draws(field$eta[column_major_cells, ])
  a <- Filter(function(d) d$scn == "A", draws)
  b <- Filter(function(d) d$scn == "B", draws)
  # Scenario A: binomial truth. Coverage 0.93, median 1,249 cells.
  res <- coverage_of(a, field, "none")
  expect_equal(res$coverage, 0.930)
  expect_equal(res$median_cells, 1249)
  expect_lt(abs(res$coverage - 0.948), 3 * sqrt(0.948 * 0.052 / 500))
  # Scenario B: a station-bin effect of SD 1. The logit-normal
  # likelihood covers 0.94 with a median of 2,491 cells; the plain binomial on
  # the population-averaged curve covers 0.658 with 1,134.
  res <- coverage_of(b, field, "station_hour", nodes = 30)
  expect_equal(res$coverage, 0.940)
  expect_equal(res$median_cells, 2491)
  gh <- gauss_hermite(30)
  averaged <- field
  averaged$eta <- stats::qlogis(Reduce(`+`, lapply(seq_along(gh$z), function(q) {
    gh$w[q] * stats::plogis(field$eta + sqrt(2) * gh$z[q])
  })))
  res <- coverage_of(b, averaged, "none")
  expect_equal(res$coverage, 0.658)
  expect_equal(res$median_cells, 1134)
})

test_that("20,000 scenario-A bins cover 0.9481", {
  skip_on_cran()
  field <- fixture_square_lake()
  column_major_cells <- match(column_major_to_terra(1:10000), field$cells)
  p <- stats::plogis(field$eta[column_major_cells, ])
  set.seed(11)
  draws <- list()
  skipped <- 0
  while (length(draws) < 20000) {
    tc <- sample.int(10000, 1)
    y <- stats::rbinom(5, 20, p[tc, ])
    if (sum(y) == 0) {
      skipped <- skipped + 1
      next
    }
    draws[[length(draws) + 1]] <- list(tc = tc, y = y)
  }
  expect_equal(skipped, 1332)
  chunks <- split(seq_along(draws), rep(1:10, each = 2000))
  hits <- vapply(chunks, function(i) coverage_of(draws[i], field, "none")$coverage * 2000,
                 numeric(1))
  coverage <- sum(hits) / 20000
  expect_equal(coverage, 0.9481)
  expect_lt(abs(coverage - 0.948), 0.0048)
})

test_that("a sigma below the true per-bin spread gives regions that miss the truth", {
  # A sigma below the true per-bin spread: 1.362, the station-hour
  # component alone, against the whole per-bin spread of 2.55.
  # Truth: logit-normal per receiver-bin, SD 2.55, n = 30.
  skip_on_cran()
  field <- fixture_square_lake()
  set.seed(20261029)
  draws <- list()
  skipped <- 0
  while (length(draws) < 1000) {
    tc <- sample.int(nrow(field$eta), 1)
    y <- stats::rbinom(5, 30, stats::plogis(field$eta[tc, ] + stats::rnorm(5, 0, 2.55)))
    if (sum(y) == 0) {
      skipped <- skipped + 1
      next
    }
    draws[[length(draws) + 1]] <- list(tc = tc, y = y)
  }
  expect_equal(skipped, 28)
  counts <- do.call(rbind, lapply(seq_along(draws), function(i) {
    counts_for(field, draws[[i]]$y, 30, "rep", i)
  }))
  start <- as.POSIXct("2024-06-01", tz = "UTC") + seq_along(draws) * 3600
  truth <- paste(start, field$cells[vapply(draws, `[[`, numeric(1), "tc")])
  coverage <- function(sigma) {
    field$model$sigma$station_hour <- sigma
    s <- location_surface(counts, field, control = list(kept_mass = 0.95))
    mean(truth %in% paste(s$bin_start, s$cell))
  }
  low <- coverage(1.362)
  right <- coverage(2.55)
  expect_equal(c(low, right), c(0.685, 0.952))
  expect_lt(low, 0.95 - 10 * sqrt(0.95 * 0.05 / 1000))
  expect_lt(abs(right - 0.95), 3 * sqrt(0.95 * 0.05 / 1000))
})

# ---- Coverage when the fitted form is wrong ----
# Every coverage test above simulates the form it fits, so none could
# see a wrong shape of overdispersion. Truth here is a hurdle: a receiver
# hears the tag in the hour with probability q(eta), and given that, logs a
# zero-truncated logit-normal count, with values from a range-window fit.
# The array is a 5 x 5 grid at 600 m, dense enough for silent receivers to
# matter.
dense_array <- function() {
  g <- terra::rast(nrows = 100, ncols = 100, xmin = 0, xmax = 3000, ymin = 0, ymax = 3000,
                   crs = "EPSG:32617", vals = 1)
  xy <- expand.grid(x = seq(315, 2715, by = 600), y = seq(315, 2715, by = 600))
  rx <- data.frame(receiver_serial = sprintf("G%02d", seq_len(nrow(xy))), x = xy$x, y = xy$y)
  fixture_field(g, rx, range_curve_eta, sigma = 2.55)
}

r_hurdle <- function(eta, n) {
  vapply(eta, function(e) {
    if (stats::runif(1) > stats::plogis(2.381 + 1.313 * e)) return(0)
    repeat {
      y <- stats::rbinom(1, n, stats::plogis(1.675 + 0.841 * e + stats::rnorm(1, 0, 1.71)))
      if (y > 0) return(y)
    }
  }, numeric(1))
}

r_logitnormal <- function(eta, n) {
  stats::rbinom(length(eta), n, stats::plogis(eta + stats::rnorm(length(eta), 0, 2.55)))
}

coverage_under <- function(field, draw_y, k = 300, n = 30) {
  set.seed(20261030)
  draws <- list()
  while (length(draws) < k) {
    tc <- sample.int(nrow(field$eta), 1)
    y <- draw_y(field$eta[tc, ], n)
    if (sum(y) > 0) draws[[length(draws) + 1]] <- list(tc = tc, y = y)
  }
  counts <- do.call(rbind, lapply(seq_along(draws), function(i) {
    counts_for(field, draws[[i]]$y, n, "rep", i)
  }))
  s <- location_surface(counts, field, control = list(kept_mass = 0.95))
  start <- as.POSIXct("2024-06-01", tz = "UTC") + seq_along(draws) * 3600
  truth <- paste(start, field$cells[vapply(draws, `[[`, numeric(1), "tc")])
  mean(truth %in% paste(s$bin_start, s$cell))
}

test_that("the station-hour form under-covers when the truth is a hurdle", {
  skip_on_cran()
  field <- dense_array()
  se <- sqrt(0.95 * 0.05 / 300)
  right <- coverage_under(field, r_logitnormal)
  wrong <- coverage_under(field, r_hurdle)
  # 280 and 191 of 300 detected bins.
  expect_equal(round(c(right, wrong), 4), c(0.9333, 0.6367))
  expect_lt(abs(right - 0.95), 3 * se)
  expect_lt(wrong, 0.95 - 10 * se)
})

# The default nodes against a 400-node, step-0.002 reference built
# through the internals, at the surface level, for n of 30, 80 and 1,440 and
# counts from 0 to n, at both sigmas tested above.
reference_surface <- function(field, y, n, sigma) {
  gh <- gauss_hermite(400)
  od <- list(form = "station_hour", sigma = sigma, nodes = 400L, z = gh$z, w = gh$w,
             step = 0.002)
  prep <- ls_prepare(field$eta, od, n)
  ll <- ls_loglik(seq_along(y), y, n, prep)
  w <- exp(ll - max(ll))
  w / sum(w)
}

test_that("default surfaces match 400 nodes to a total variation below 1e-4", {
  skip_on_cran()
  rx <- data.frame(receiver_serial = c("R1", "R2", "R3", "R4"), x = c(150, 450, 150, 450),
                   y = c(150, 150, 450, 450))
  field <- fixture_field(fixture_grid(20), rx, function(d) 2 - d / 100)
  tv <- c()
  for (sigma in c(2.55, 3.4)) {
    for (n in c(30, 80, 1440)) {
      for (y in list(c(n, round(n / 2), 1, 0), c(0, 0, 3, 0), c(round(n / 3), n, 0, 7))) {
        field$model$sigma$station_hour <- sigma
        got <- full_surface(counts_for(field, y, n), field)
        ref <- reference_surface(field, y, n, sigma)
        tv <- c(tv, 0.5 * sum(abs(got - ref)))
        expect_equal(length(region95(got)), length(region95(ref)))
      }
    }
  }
  expect_length(tv, 18)
  expect_lt(max(tv), 1e-4)
})

test_that("the cubic interpolation between table nodes matches direct quadrature", {
  skip_on_cran()
  set.seed(4)
  for (sigma in c(2.55, 3.4)) {
    od <- ls_overdispersion("station_hour", sigma, ls_control(list()),
                            data.table::data.table(bin_start = 0, bin_end = 3600))
    for (n in c(30, 1440)) {
      # Logits between the nodes of the step-0.01 table, none on a node.
      eta <- sort(stats::runif(200, -12, 6))
      eta <- eta[abs(eta / 0.01 - round(eta / 0.01)) > 0.1]
      prep <- ls_prepare(matrix(eta), od, n)
      for (y in unique(c(0, 1, round(n / 2), n))) {
        direct <- ls_table(c(eta, max(eta) + 1), y, n, od)$c0
        expect_lt(max(abs(ls_term(eta, y, n, prep) - direct)), 1e-6)
      }
    }
  }
})

# nolint end
