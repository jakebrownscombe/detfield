#' A logistic detection curve from two distances
#'
#' Builds a [de_model()] whose detection probability falls with distance
#' along a logistic curve that passes through 0.5 at `d50` and 0.05 at `d05`.
#' On the logit scale the curve is `qlogis(0.05) * (distance - d50) / (d05 -
#' d50)`, so the probability at the receiver is above 0.5 whenever `d50` is
#' above 0.
#'
#' With `depth`, each of `d50` and `d05` holds one value per receiver depth.
#' Between two depths both distances are interpolated linearly; outside the
#' range they are held at the nearest end.
#'
#' positionR 1.x called the 5 per cent distance `d95`. It is renamed `d05`
#' here because the curve reaches 0.05, not 0.95, there.
#'
#' @param d50 Distance in metres at which half of the transmissions are
#'   detected; at least 0.
#' @param d05 Distance in metres at which 5 per cent are detected; larger
#'   than `d50`.
#' @param depth Optional receiver depths in metres, increasing, one per value
#'   of `d50` and `d05`.
#'
#' @return A `de_model` with no `sigma` and no modulators. Use [de_model()] on
#'   it to add them.
#'
#' @seealso [de_model()], [de_field()].
#' @export
#' @examples
#' # The Stoney curve at a 1 m receiver: half heard at 250 m, 5% at 500 m.
#' m <- de_logistic(d50 = 250, d05 = 500)
#' plogis(m$predict(c(0, 250, 500)))
#'
#' # The Stoney example: a longer range at a deeper receiver, linear in depth.
#' m_depth <- de_logistic(d50 = c(250, 750), d05 = c(500, 1500), depth = c(1, 35))
#' plogis(m_depth$predict(c(500, 500), depth_m = c(1, 18)))
#'
#' # Add the per-bin spread, estimated around this same curve.
#' de_model(m, sigma = list(station_hour = 1))
de_logistic <- function(d50, d05, depth = NULL) {
  depth <- check_logistic(d50, d05, depth)
  by_depth <- length(depth) > 1
  # Four significant digits: a fitted d50 printed to 15 digits is unreadable.
  short <- function(v) paste(signif(v, 4), collapse = "/")
  src <- if (by_depth) {
    sprintf("de_logistic, d50 %s m, d05 %s m, depths %s m", short(d50), short(d05),
            short(depth))
  } else {
    sprintf("de_logistic, d50 %s m, d05 %s m", short(d50), short(d05))
  }
  core <- list(predict = logistic_curve(d50, d05, depth), source = src,
               uses = c("distance_m", if (by_depth) "depth_m"))
  new_de_model(core, check_sigma(list()), list())
}

# The depths the curve is tabulated at; 0 stands for "no depth".
check_logistic <- function(d50, d05, depth) {
  distances_ok <- function(v) is.numeric(v) && length(v) && all(is.finite(v)) && all(v >= 0)
  if (!distances_ok(d50) || !distances_ok(d05)) {
    stop("`d50` and `d05` must be finite numbers of at least 0.", call. = FALSE)
  }
  if (length(d50) != length(d05)) stop("`d50` and `d05` must have the same length.", call. = FALSE)
  if (any(d05 <= d50)) {
    stop("`d05` must be larger than `d50`: detection falls from 0.5 to 0.05 with distance.",
         call. = FALSE)
  }
  check_logistic_depth(depth, length(d50))
}

check_logistic_depth <- function(depth, n) {
  if (is.null(depth)) {
    if (n != 1) stop("Give `depth` when `d50` has more than one value.", call. = FALSE)
    return(0)
  }
  ok <- is.numeric(depth) && length(depth) == n && all(is.finite(depth)) && all(diff(depth) > 0)
  if (!ok) stop("`depth` must be increasing finite numbers, one per value of `d50`.", call. = FALSE)
  depth
}

logistic_curve <- function(d50, d05, depth) {
  by_depth <- length(depth) > 1
  function(distance_m, depth_m = NA_real_, tag_model = NA_character_, ...) {
    if (!is.numeric(distance_m)) stop("`distance_m` must be numeric.", call. = FALSE)
    z <- rep_len(as.numeric(depth_m), length(distance_m))
    if (by_depth && anyNA(z)) stop("This curve needs `depth_m`.", call. = FALSE)
    a <- if (by_depth) stats::approx(depth, d50, xout = z, rule = 2)$y else d50
    b <- if (by_depth) stats::approx(depth, d05, xout = z, rule = 2)$y else d05
    stats::qlogis(0.05) * (distance_m - a) / (b - a)
  }
}
