# Title: Error of one receiver's term under the station-hour spread
# Description: The numbers that ?location_surface states under Computation
#   about the error of each receiver's term, and the comment above
#   ls_nodes() in R/likelihood.R. One receiver's term,
#   log E[Binomial(y; n, plogis(eta + sigma * Z))] with Z standard normal,
#   is computed the way location_surface() computes it, from the table and
#   its cubic interpolation, and compared with stats::integrate(). The table
#   starts below the smallest logit here, so the closed form that
#   location_surface() uses below the table is not tested by this script.
#   Counts y run from 0 to n, for n of 30, 60 and 720 transmissions, at 39
#   logits from -19.995 to 7.745. Each logit sits half a table step from the
#   nearest grid point, where the interpolation is farthest from a table value.
# Run from the package root: Rscript validation/02_quadrature_error.R
# Output: validation/outputs/02_quadrature_error.csv and the printed lines
#   below, kept in validation/outputs/02_quadrature_error.txt.

# ---- Libraries ----
library(detfield)

# ---- Settings ----
sigmas <- c(2.55, 3.4, 6)
n_values <- c(30, 60, 720)
step <- 0.01  # the default control$table_step
# The table's grid points are whole multiples of the step, so a shift of half
# a step puts every logit midway between two of them.
eta <- seq(-20, 8, by = 0.73) + step / 2
length(eta)  # logits per count
range(abs(eta / step - round(eta / step)))  # distance from the grid, in steps

# ---- The term as location_surface() computes it ----
term_detfield <- function(eta, y, n, sigma, nodes) {
  gh <- detfield:::gauss_hermite(nodes)
  od <- list(form = "station_hour", sigma = sigma, nodes = nodes, z = gh$z, w = gh$w,
             step = step)
  prep <- detfield:::ls_prepare(matrix(eta), od, n)
  detfield:::ls_term(eta, y, n, prep)
}

# ---- The term by integrate() ----
# The integrand in z is log-concave, so it has one mode. integrate() runs on
# each side of the mode, out to 40 standard normal units, where the normal
# factor alone is below exp(-800) of its peak.
term_exact <- function(eta, y, n, sigma) {
  vapply(eta, function(e) {
    # The binomial log probability in the logit, so it stays finite where
    # plogis() rounds to 0 or 1.
    f <- function(z) {
      u <- e + sigma * z
      lchoose(n, y) + y * stats::plogis(u, log.p = TRUE) +
        (n - y) * stats::plogis(-u, log.p = TRUE) + stats::dnorm(z, log = TRUE)
    }
    m <- stats::optimize(f, c(sigma * (y - n) - 1, sigma * y + 1), maximum = TRUE,
                         tol = 1e-10)
    top <- m$objective
    g <- function(z) exp(f(z) - top)
    left <- stats::integrate(g, m$maximum - 40, m$maximum, rel.tol = 1e-12,
                             subdivisions = 2000L)$value
    right <- stats::integrate(g, m$maximum, m$maximum + 40, rel.tol = 1e-12,
                              subdivisions = 2000L)$value
    top + log(left + right)
  }, numeric(1))
}

# ---- Error over counts and logits ----
# Nodes: the default rule (25 per logit of sigma, 20 to 200) and 20 nodes.
res <- do.call(rbind, lapply(sigmas, function(sigma) {
  rule <- detfield:::ls_nodes(sigma)
  do.call(rbind, lapply(n_values, function(n) {
    err <- vapply(c(rule, 20), function(nodes) {
      max(vapply(0:n, function(y) {
        exact <- term_exact(eta, y, n, sigma)
        ok <- is.finite(exact)
        max(abs(term_detfield(eta, y, n, sigma, nodes)[ok] - exact[ok]))
      }, numeric(1)))
    }, numeric(1))
    data.frame(sigma = sigma, n = n, rule_nodes = rule, max_err_rule = err[1],
               max_err_20 = err[2])
  }))
}))
print(res, digits = 3)

# ---- Worst case per standard deviation ----
# Over the three numbers of transmissions; each row covers 39 logits and
# every count from 0 to n.
worst <- aggregate(cbind(max_err_rule, max_err_20) ~ sigma + rule_nodes, data = res, FUN = max)
print(worst, digits = 3)
write.csv(res, "validation/outputs/02_quadrature_error.csv", row.names = FALSE)

# ---- Session ----
sessionInfo()
