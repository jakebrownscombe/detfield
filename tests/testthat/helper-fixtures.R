# ---- Fixtures for the likelihood tests ----
# A detection field in the layout location_surface() reads (eta, cells, receivers,
# grid, model), built here so these tests do not depend on de_field().

# A square grid of `side` by `side` cells, all water unless `land` lists cells.
fixture_grid <- function(side = 10, cell = 30, land = integer()) {
  v <- rep(1, side * side)
  v[land] <- NA
  terra::rast(nrows = side, ncols = side, xmin = 0, xmax = side * cell, ymin = 0,
              ymax = side * cell, crs = "EPSG:32617", vals = v)
}

# eta_fun maps a matrix of distances in metres to logits. The default sigma,
# 2.55 logits, is a per-hour spread of the size a lake range test gives.
fixture_field <- function(grid, rx, eta_fun, sigma = 2.55, key = NULL) {
  cells <- which(!is.na(terra::values(grid, mat = FALSE)))
  xy <- terra::xyFromCell(grid, cells)
  d <- sqrt(outer(xy[, 1], rx$x, "-")^2 + outer(xy[, 2], rx$y, "-")^2)
  ids <- rx[[if (is.null(key)) "receiver_serial" else key]]
  structure(list(
    eta = eta_fun(d), cells = cells, receivers = ids, grid = grid, receiver_xy = rx,
    key = key,
    model = structure(list(
      predict = function(...) stop("the likelihood must not call the model"),
      sigma = list(station_hour = sigma, station_day = NA), modulators = list(),
      source = "test fixture"
    ), class = "de_model")
  ), class = "de_field")
}

fixture_receivers <- function() {
  data.frame(receiver_serial = c("R1", "R2", "R3", "R4"),
             station_id = c("S1", "S2", "S3", "S4"),
             x = c(75, 225, 75, 225), y = c(75, 75, 225, 225))
}

# A logistic range curve at 10 m station depth, with the full-precision
# coefficients of a quasibinomial range model on distance and station depth.
range_curve_eta <- function(d) {
  -0.33820268214301519 - 0.0032194253046818667 * d + 0.10594944808292814 * 10
}

# A square known-truth lake: 100 by 100 cells of 30 m, five stations by default.
fixture_square_lake <- function(stations = 5, eta_fun = range_curve_eta, sigma = 1) {
  grid <- terra::rast(nrows = 100, ncols = 100, xmin = 0, xmax = 3000, ymin = 0,
                      ymax = 3000, crs = "EPSG:32617", vals = 1)
  rx <- data.frame(receiver_serial = c("C", "N", "S", "E", "W", "NE", "NW", "SE", "SW"),
                   x = c(1515, 1515, 1515, 2115, 915, 2715, 315, 2715, 315),
                   y = c(1515, 2115, 915, 1515, 1515, 2715, 2715, 315, 315))
  fixture_field(grid, rx[seq_len(stations), ], eta_fun, sigma = sigma)
}

# A grid built by data.table::CJ(col, row) with cell = .I numbers its cells
# column by column: cell k sits at column ceiling(k / 100) and row
# k - 100 * (column - 1), rows counted up from the bottom. This maps such a
# cell number to a terra cell.
column_major_to_terra <- function(k) {
  col <- (k - 1) %/% 100 + 1
  row <- (k - 1) %% 100 + 1
  (100 - row) * 100 + col
}

counts_for <- function(field, y, n, animal = "a1", start = 0, key = "receiver_serial") {
  t0 <- as.POSIXct("2024-06-01 00:00:00", tz = "UTC") + start * 3600
  out <- data.frame(animal_id = animal, bin_start = t0, bin_end = t0 + 3600,
                    receiver_serial = field$receivers, y = as.integer(y), n = as.integer(n),
                    listening = TRUE, stringsAsFactors = FALSE)
  if (key == "station_id") out$station_id <- field$receivers
  out
}

# Brute-force Bayes: one loop over cells and one over receivers, dbinom as R
# computes it, a flat prior, and normalisation over water cells.
brute_posterior <- function(field, y, n) {
  p <- stats::plogis(field$eta)
  ll <- numeric(nrow(p))
  for (c in seq_len(nrow(p))) {
    for (s in seq_len(ncol(p))) ll[c] <- ll[c] + stats::dbinom(y[s], n, p[c, s], log = TRUE)
  }
  w <- exp(ll - max(ll))
  list(post = w / sum(w), log_lik = max(ll) + log(sum(w)) - log(length(ll)))
}

# A plain logit-normal likelihood, written apart from the package: fixed
# Gauss-Hermite nodes on the probability scale, the node count an argument.
fixed_node_loglik_logitnormal <- function(y, n, eta, sigma_u, nodes = 30) {
  k <- 1:(nodes - 1)
  jac <- matrix(0, nodes, nodes)
  jac[cbind(k, k + 1)] <- sqrt(k / 2)
  jac[cbind(k + 1, k)] <- sqrt(k / 2)
  e <- eigen(jac, symmetric = TRUE)
  gh <- list(z = e$values, w = e$vectors[1, ]^2)
  ll <- numeric(nrow(eta))
  for (s in seq_along(y)) {
    acc <- 0
    for (q in seq_along(gh$z)) {
      acc <- acc + gh$w[q] * stats::dbinom(y[s], n, stats::plogis(eta[, s] + sqrt(2) * sigma_u *
                                                                    gh$z[q]))
    }
    ll <- ll + log(pmax(acc, 1e-300))
  }
  ll
}

# The full posterior of one bin, kept_mass = 1, ordered by cell.
full_surface <- function(counts, field, form = "station_hour", ...) {
  s <- suppressMessages(location_surface(counts, field, form,
                                         control = list(kept_mass = 1, ...)))
  v <- numeric(max(field$cells))
  v[s$cell] <- s$prob
  v[field$cells]
}

region95 <- function(prob) {
  o <- order(prob, decreasing = TRUE)
  o[seq_len(which(cumsum(prob[o]) >= 0.95)[1])]
}
