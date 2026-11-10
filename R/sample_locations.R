#' Draw locations from location surfaces
#'
#' @description
#' Draws cells from location surfaces in proportion to their probability,
#' with replacement, for analyses that take points rather than surfaces, such
#' as habitat selection.
#'
#' @details
#' Draws are with replacement, so a cell is drawn as often as its probability
#' says; a draw without replacement would push draws out of the most probable
#' cells once those were taken. Each draw is the centre of its cell. Only the
#' kept cells of a surface can be drawn (see [location_surface()]); their
#' probabilities are rescaled to sum to one. `n` cells are drawn for each
#' animal and bin, so the draws describe the animal's location in its
#' detected bins only.
#'
#' @param surface A location surface from [location_surface()].
#' @param n Number of draws per animal and bin.
#' @param seed An integer seed for the draws, or `NULL` to use the current
#'   random number stream. The caller's stream is restored afterwards when a
#'   seed is given.
#' @return A data frame with one row per draw: `animal_id`, `bin_start` (the
#'   bin the draw came from), `draw`, `cell`, `x` and `y` (the cell centre in
#'   the grid's coordinates, `NA` when the surface holds no grid). Attributes:
#'   `"estimand"` and `"epsg"`, the grid's EPSG code.
#' @seealso [location_surface()], [credible_region()]
#' @export
#' @examples
#' tif <- system.file("extdata", "stoney_depth.tif", package = "detfield")
#' grid <- terra::rast(tif)
#' # Beyond 2 km the Stoney curve is below 0.001 at every receiver; leaving
#' # those pairs out keeps the example quick.
#' dist <- water_distances(grid, stoney_track$receivers,
#'                         max_distance = 2000)
#' model <- de_logistic(d50 = c(250, 750), d05 = c(500, 1500), depth = c(1, 35))
#' field <- de_field(dist, de_model(model, sigma = list(station_hour = 0)))
#' counts <- bin_detections(stoney_track$detections, stoney_track$receivers,
#'                          stoney_track$tags)
#' surface <- location_surface(counts, field)
#' head(sample_locations(surface, n = 100, seed = 1))
sample_locations <- function(surface, n, seed = NULL) {
  sl_check_surface(surface)
  if (!ls_whole(n) || n < 1) stop("`n` must be one whole number of at least 1.", call. = FALSE)
  grid <- sl_grid(surface, strict = TRUE)
  s <- data.table::data.table(animal_id = surface$animal_id, bin_start = surface$bin_start,
                              cell = surface$cell, prob = surface$prob)
  # The result is built under the seed too, so a seeded call leaves the
  # session's random stream as it found it.
  with_seed_or_not(seed, sl_result(sl_draw(s, as.integer(n)), grid))
}

# Draws from every bin in one pass: the weights are normalised within bins,
# bin g occupies the interval (g - 1, g] of one cumulative scale, and a
# uniform draw in (g - 1, g) falls in exactly one of its rows.
sl_draw <- function(s, n) {
  prob <- w <- g <- draw <- NULL
  group <- c("animal_id", "bin_start")
  s[, w := prob / sum(prob), by = group]
  data.table::setorderv(s, c(group, "cell"))
  s[, g := data.table::rleidv(s, group)]
  edge <- s$g - 1 + s[, cumsum(w) / sum(w), by = "g"]$V1
  edge[!duplicated(s$g, fromLast = TRUE)] <- s$g[!duplicated(s$g, fromLast = TRUE)]
  k <- max(s$g)
  u <- rep(seq_len(k) - 1, each = n) + stats::runif(k * n)
  rows <- findInterval(u, edge) + 1L
  out <- s[rows, c("animal_id", "bin_start", "cell")]
  out[, draw := rep(seq_len(n), times = k)]
  out
}

sl_result <- function(drawn, grid) {
  xy <- if (is.null(grid)) {
    matrix(NA_real_, nrow(drawn), 2)
  } else {
    terra::xyFromCell(grid, drawn$cell)
  }
  out <- data.frame(animal_id = drawn$animal_id, bin_start = drawn$bin_start,
                    draw = drawn$draw, cell = drawn$cell, x = xy[, 1], y = xy[, 2],
                    stringsAsFactors = FALSE)
  attr(out, "estimand") <- "location of the animal in each detected bin"
  attr(out, "epsg") <- grid_epsg(grid)
  out
}

# ---- Helpers shared with credible_region() and the methods ----
sl_check_surface <- function(surface) {
  if (!inherits(surface, "loc_surface")) {
    stop("`surface` must be a location surface from location_surface().", call. = FALSE)
  }
  if (!nrow(surface)) stop("`surface` holds no cells.", call. = FALSE)
  invisible(surface)
}

# The surface's grid as a SpatRaster, or NULL when it holds none. A grid
# stored as a path that cannot be read gives a warning, or with strict =
# TRUE an error, rather than NA coordinates without notice.
sl_grid <- function(surface, strict = FALSE) {
  g <- attr(surface, "grid")
  if (is.null(g) || inherits(g, "SpatRaster")) return(g)
  if (inherits(g, "PackedSpatRaster")) return(terra::rast(g))
  readable <- is.character(g) && length(g) == 1 && file.exists(g)
  if (readable) return(terra::rast(g))
  msg <- sprintf("The surface's grid (%s) cannot be read from %s.",
                 if (is.character(g)) g[1] else class(g)[1], getwd())
  if (strict) stop(msg, call. = FALSE)
  warning(paste(msg, "Cell areas and coordinates are NA."), call. = FALSE)
  NULL
}

# Evaluates `code` under `seed` and leaves the session's random stream as it
# found it; with no seed, `code` draws from the session's stream.
with_seed_or_not <- function(seed, code) {
  if (is.null(seed)) code else withr::with_seed(seed, code)
}
